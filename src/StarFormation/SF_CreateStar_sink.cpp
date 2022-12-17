#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_AGORA
// Description :  Create new star particles stochastically using the presription suggested by the AGORA project
//
// Note        :  1. Ref: (1) Nathan Goldbaum, et al., 2015, ApJ, 814, 131 (arXiv: 1510.08458), sec. 2.4
//                        (2) Ji-hoon Kim, et al., 2016, ApJ, 833, 202 (arXiv: 1610.03066), sec. 3.2
//                2. One must turn on STORE_POT_GHOST when adopting STORE_PAR_ACC
//                   --> It is because, currently, this function always uses the pot_ext[] array of each patch
//                       to calculate the gravitationally acceleration of the new star particles
//                3. One must invoke Buf_GetBufferData( ..., _TOTAL, ... ) after calling this function
//                4. Currently this function does not check whether the cell mass exceeds the Jeans mass
//                   --> Ref: "jeanmass" in star_maker_ssn.F of Enzo
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Current physical time (after advancing solution by dt)
//                dt           : Time interval to advance solution
//                               --> Currently this function does not distinguish dt and the physical time interval (dTime)
//                               --> Does NOT support COMOVING yet
//                RNG          : Random number generator
//                GasDensThres : Minimum gas density for creating star particles                (--> "SF_CREATE_STAR_MIN_GAS_DENS"  )
//                Efficiency   : Gas-to-star mass efficiency                                    (--> "SF_CREATE_STAR_MASS_EFF"      )
//                MinStarMass  : Minimum star particle mass for the stochastical star formation (--> "SF_CREATE_STAR_MIN_STAR_MASS" )
//                MaxStarMFrac : Maximum gas mass fraction allowed to convert to stars          (--> "SF_CREATE_STAR_MAX_STAR_MFRAC")
//                DetRandom    : Make random numbers determinisitic                             (--> "SF_CREATE_STAR_DET_RANDOM"    )
//                UseMetal     : Store the metal mass fraction in star particles
//
// Return      :  1. Particle repository will be updated
//                2. fluid[] array of gas will be updated
//-------------------------------------------------------------------------------------------------------
void SF_CreateStar_AGORA( const int lv, const real TimeNew, const real dt, RandomNumber_t *RNG,
                          const real GasDensThres, const real Efficiency, const real MinStarMass, const real MaxStarMFrac,
                          const bool DetRandom, const bool UseMetal, const real AccCellNum)
{

// check
#  if ( defined STORE_PAR_ACC  &&  !defined STORE_POT_GHOST )
#     error : STAR_FORMATION + STORE_PAR_ACC must work with STORE_POT_GHOST !!
#  endif

#  ifndef GRAVITY
#     error : must turn on GRAVITY for SF_CreateStar_sink() !!
#  endif

#  ifdef COMOVING
#     error : SF_CreateStar_sink() does not support COMOVING yet !!
#  endif

#  ifdef GAMER_DEBUG
   if ( UseMetal  &&  Idx_ParMetalFrac == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParMetalFrac is undefined for \"UseMetal\" !!\n" );

   if ( UseMetal  &&  Idx_Metal == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_Metal is undefined for \"UseMetal\" !!\n" );

   if ( Idx_ParCreTime == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParCreTime is undefined !!\n" );
#  endif // #ifdef GAMER_DEBUG

// checking the value of accretion radius
   if (AccCellNum > PS1)
      Aux_Error( ERROR_INFO, "AccCellNum should be smaller than PATCH_SIZE !!" );

// constant parameters
   const double dh             = amr->dh[lv];
   
   const double AccRadius      = AccCellNum*dh;  /////////////////////////////////
   const int    NGhost         = AccCellNum;
   const int    Size_Flu       = PS2 + 2*NGhost; // final cube size
   const int    Size_Flu_P1    = Size_Flu + 1; // for face-centered B field
   const int    Size_Pot       = Size_Flu; // for potential
   const int    NPG            = 1;

   const real   dv             = CUBE( dh );
   const int    FluSg          = amr->FluSg[lv];
   const int    PotSg          = amr->PotSg[lv];
   const real   Coeff_FreeFall = SQRT( (32.0*NEWTON_G)/(3.0*M_PI) );
   const real  _MinStarMass    = (real)1.0 / MinStarMass;
   const real   Eff_times_dt   = Efficiency*dt;
// const real   GraConst       = ( OPT__GRA_P5_GRADIENT ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh);
   const real   GraConst       = ( false                ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh); // P5 is NOT supported yet


// start of OpenMP parallel region
#  pragma omp parallel
   {

// thread-private variables
#  ifdef OPENMP
   const int TID = omp_get_thread_num();
#  else
   const int TID = 0;
#  endif

   double x0, y0, z0, x, y, z, VelX, VelY, VelZ;
   real   GasDens, GasDensFreeFall, _GasDens, GasMass, _Time_FreeFall, StarMFrac, StarMass, GasMFracLeft;
   real   (*fluid)[PS1][PS1][PS1]      = NULL; // fluid pointer, PS1: PATCH_SIZE
#  ifdef STORE_POT_GHOST
   real   (*pot_ext)[GRA_NXT][GRA_NXT] = NULL;
#  endif

   const int MaxNewParPerPatch = CUBE(PS1);
   real   (*NewParAtt)[PAR_NATT_TOTAL] = new real [MaxNewParPerPatch][PAR_NATT_TOTAL];
   long    *NewParID                   = new long [MaxNewParPerPatch];

   int NNewPar;

// load the existing particles ID (the number)
   int   NParTot   = amr->Par->NPar_Active_AllRank;
   const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const real *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

// loop over all real patch groups
// use static schedule to ensure bitwise reproducibility when running with the same numbers of OpenMP threads and MPI ranks
// --> bitwise reproducibility will still break when running with different numbers of OpenMP threads and/or MPI ranks
//     unless both BITWISE_REPRODUCIBILITY and SF_CREATE_STAR_DET_RANDOM are enabled
#  pragma omp for schedule( static )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
      real (*Flu_Array_F_In)[CUBE(Size_Flu)] = new real [FLU_NIN][CUBE(Size_Flu)];
      real (*Mag_Array_F_In)[Size_Flu_P1*SQR(Size_Flu)] = new real [NCOMP_MAG][Size_Flu_P1*SQR(Size_Flu)];
      real (*Pot_Array_USG_F) = new real [CUBE(Size_Pot)];
      real fluid[FLU_NIN];
      real Corner_Array_F[3]; // the corner of the ghost zone

#     if ( MODEL != HYDRO )
      const double MIN_DENS            = -1.0;  // set to an arbitrarily negative value to disable it
#     endif
#     ifndef MHD
      const int    OPT__MAG_INT_SCHEME = INT_NONE;
#     endif
      const bool   IntPhase_No         = false;
      const real   MinDens_No          = -1.0;
      const real   MinPres_No          = -1.0;
      const real   MinTemp_No          = -1.0;
      const real   MinEntr_No          = -1.0;
      const bool   DE_Consistency_Yes  = true;
      const bool   DE_Consistency_No   = false;
      const bool   DE_Consistency      = ( OPT__OPTIMIZE_AGGRESSIVE ) ? DE_Consistency_No : DE_Consistency_Yes;
      const real   MinDens             = ( OPT__OPTIMIZE_AGGRESSIVE ) ? MinDens_No : MIN_DENS;

#     ifdef MHD
      real *Mag_Array = Mag_Array_F_In[0][0];
#     else
      real *Mag_Array = NULL;
#     endif
      Prepare_PatchData( lv, TimeNew, Flu_Array_F_In[0][0], Mag_Array,
                        NGhost, NPG, &PID0, _TOTAL, _MAG,
                        OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                        OPT__BC_FLU, BC_POT_NONE, MinDens,    MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency );

#     ifdef UNSPLIT_GRAVITY
//    prepare the potential array
      if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
      Prepare_PatchData( lv, TimeNew, Pot_Array_USG_F[0], NULL,
                        NGhost, NPG, &PID0, _POTE, _NONE,
                        OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                        OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

//    prepare the corner array
     if ( OPT__EXT_ACC )
      {
         for (int d=0; d<3; d++)    Corner_Array_F[d] = amr->patch[0][lv][PID0]->EdgeL[d] + 0.5*dh - dh*NGhost;
      }
#     endif // #ifdef UNSPLIT_GRAVITY

      for (int pk=NGhost; pk<PS2 + NGhost; k++)
      for (int pj=NGhost; pj<PS2 + NGhost; j++)
      for (int pi=NGhost; pi<PS2 + NGhost; i++) // loop inside the patch group
      {
         x = Corner_Array_F[0] + pi*dh + dh*NGhost;
         y = Corner_Array_F[1] + pj*dh + dh*NGhost;
         z = Corner_Array_F[2] + pk*dh + dh*NGhost;

         const int t = IDX321( pi, pj, pk, FLU_NXT, FLU_NXT );
         for (int v=0; v<FLU_NIN; v++)    fluid[v] = Flu_Array_F_In[v][t];
         VelX = fluid[MOMX]/fluid[DENS];
         VelY = fluid[MOMY]/fluid[DENS];
         VelZ = fluid[MOMZ]/fluid[DENS];

//       1. Proximity Check + Density Threshold
//       ===========================================================================================================
         GasDens = fluid[DENS];
         real GasDensThresCopy = GasDensThres;
         bool InsideAccRadius = false;

         for (int p=0; p<NParTot; p++)
         {
            real PCPX, PCPY, PCPZ; // particle-cell relative position
            real PCVX, PCVY, PCVZ; // particle-cell relative velocity
            real D2Par; // particle-cell distance

            PCPX = x - ParPos[0][p];
            PCPY = y - ParPos[1][p];
            PCPZ = z - ParPos[2][p];

            D2Par = SQRT(SQR(PCPX)+SQR(PCPX)+SQR(PCPX));
            if ( D2Par < AccRadius )
            {
               InsideAccRadius = true;
               break;
            }

            PCVX = VelX - ParVel[0][p];
            PCVY = VelY - ParVel[1][p];
            PCVZ = VelZ - ParVel[2][p];

            real NPCPX, NPCPY, NPCPZ; // normalized particle-cell relative position
            NPCPX = PCPX/D2Par;
            NPCPY = PCPY/D2Par;
            NPCPZ = PCPZ/D2Par;

            GasDensFreeFall = SQR((1/Coeff_FreeFall)*(NPCPX*PCVX + NPCPY*PCVY + NPCPZ*PCVZ)/D2Par); // Clarke et al. 2017, eqn (5)
            GasDensThresCopy = MAX(GasDensThresCopy, GasDensFreeFall);
         } // NParTot

         if ( InsideAccRadius )               continue;
         if ( GasDens < GasDensThresCopy )    continue;
         


         

      } // i, j, k

   }













   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip non-leaf patches
      if ( amr->patch[0][lv][PID]->son != -1 )  continue;


//    to get deterministic and different random numbers for all patches, reset the random seed of each patch according to
//    its location and time
//    --> patches at different time and/or AMR levels may still have the same random seeds...
      if ( DetRandom )
      {
//       the factor "1.0e6" in the end is just to make random seeds at different times more different, especially for
//       extremely small time-step
         const long RSeed = SF_CREATE_STAR_RSEED + amr->patch[0][lv][PID]->LB_Idx + long(TimeNew*UNIT_T/Const_yr*1.0e6);
         RNG->SetSeed( TID, RSeed );
      }


      fluid   = amr->patch[FluSg][lv][PID]->fluid;
#     ifdef STORE_POT_GHOST
      pot_ext = amr->patch[PotSg][lv][PID]->pot_ext;
#     endif
      x0      = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
      y0      = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
      z0      = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;
      NNewPar = 0;

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++)
      {

//       1. Density Threshold
//       ===========================================================================================================
         GasDens = fluid[DENS][k][j][i];
         GasMass = GasDens*dv;

         if ( GasDens < GasDensThres )    continue;

//       2. Refinement Check
//       ===========================================================================================================
//       already written in SF_CreateStar.cpp (SF_CREATE_STAR_MIN_LEVEL)

//       3. Converging Flow Check
//       ===========================================================================================================













//       1. check the star formation criteria
//       ===========================================================================================================
         GasDens = fluid[DENS][k][j][i];
         GasMass = GasDens*dv;

//       1-1. create star particles only if the gas density exceeds the given threshold
         if ( GasDens < GasDensThres )    continue;


//       1-2. estimate the gas free-fall time
//       --> consider only the gas density under the assumption that the dark matter doesn't collapse
         _Time_FreeFall = Coeff_FreeFall * SQRT( GasDens );


//       1-3. estimate the gas mass fraction to convert to stars
         StarMFrac = Eff_times_dt*_Time_FreeFall;
         StarMass  = GasMass*StarMFrac;


//       1-4. stochastic star formation
//       --> if the star particle mass (StarMass) is below the minimum mass (MinStarMass), we create a
//           new star particle with a mass of MinStarMass and a probability of StarMass/MinStarMass
//       --> Eq. [5] in Goldbaum et al. (2015)
         if ( StarMass < MinStarMass )
         {
            const double Min = 0.0;
            const double Max = 1.0;

            double Random = RNG->GetValue( TID, Min, Max );

            if ( (real)Random < StarMass*_MinStarMass )  StarMFrac = MinStarMass / GasMass;
            else                                         continue;
         }


//       1-5. check the maximum gas mass fraction allowed to convert to stars
         StarMFrac = MIN( StarMFrac, MaxStarMFrac );
         StarMass  = GasMass*StarMFrac;



//       2. store the information of new star particles
//       --> we will not create these new particles until looping over all cells in a patch in order to reduce
//           the OpenMP synchronization overhead
//       ===========================================================================================================
//       check
#        ifdef GAMER_DEBUG
         if ( NNewPar >= MaxNewParPerPatch )
            Aux_Error( ERROR_INFO, "NNewPar (%d) >= MaxNewParPerPatch (%d) !!\n", NNewPar, MaxNewParPerPatch );
#        endif

//       2-1. intrinsic attributes
         _GasDens = (real)1.0 / GasDens;
         x        = x0 + i*dh;
         y        = y0 + j*dh;
         z        = z0 + k*dh;

         NewParAtt[NNewPar][PAR_MASS] = StarMass;
         NewParAtt[NNewPar][PAR_POSX] = x;
         NewParAtt[NNewPar][PAR_POSY] = y;
         NewParAtt[NNewPar][PAR_POSZ] = z;
         NewParAtt[NNewPar][PAR_VELX] = fluid[MOMX][k][j][i]*_GasDens;
         NewParAtt[NNewPar][PAR_VELY] = fluid[MOMY][k][j][i]*_GasDens;
         NewParAtt[NNewPar][PAR_VELZ] = fluid[MOMZ][k][j][i]*_GasDens;
         NewParAtt[NNewPar][PAR_TIME] = TimeNew;
         NewParAtt[NNewPar][PAR_TYPE] = PTYPE_STAR;

//       particle acceleration
#        ifdef STORE_PAR_ACC
         real GasAcc[3] = { (real)0.0, (real)0.0, (real)0.0 };

//       external acceleration
         if ( OPT__EXT_ACC )  CPUExtAcc_Ptr( GasAcc, x, y, z, TimeNew, ExtAcc_AuxArray );

//       self-gravity and external potential
         if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
         {
            const int ii = i + GRA_GHOST_SIZE;
            const int jj = j + GRA_GHOST_SIZE;
            const int kk = k + GRA_GHOST_SIZE;

#           ifdef STORE_POT_GHOST
            const real pot_xm = pot_ext[kk  ][jj  ][ii-1];
            const real pot_xp = pot_ext[kk  ][jj  ][ii+1];
            const real pot_ym = pot_ext[kk  ][jj-1][ii  ];
            const real pot_yp = pot_ext[kk  ][jj+1][ii  ];
            const real pot_zm = pot_ext[kk-1][jj  ][ii  ];
            const real pot_zp = pot_ext[kk+1][jj  ][ii  ];
#           endif

            GasAcc[0] += GraConst*( pot_xp - pot_xm );
            GasAcc[1] += GraConst*( pot_yp - pot_ym );
            GasAcc[2] += GraConst*( pot_zp - pot_zm );
         }

         NewParAtt[NNewPar][PAR_ACCX] = GasAcc[0];
         NewParAtt[NNewPar][PAR_ACCY] = GasAcc[1];
         NewParAtt[NNewPar][PAR_ACCZ] = GasAcc[2];
#        endif // ifdef STORE_PAR_ACC


//       2-2. extrinsic attributes
//       note that we store the metal mass **fraction** instead of density in particles
         if ( UseMetal )
         NewParAtt[NNewPar][Idx_ParMetalFrac] = fluid[Idx_Metal][k][j][i] * _GasDens;

         NewParAtt[NNewPar][Idx_ParCreTime  ] = TimeNew;

         NNewPar ++;



//       3. remove the gas that has been converted to stars
//       ===========================================================================================================
         GasMFracLeft = (real)1.0 - StarMFrac;

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v][k][j][i] *= GasMFracLeft;
      } // i,j,k



//    4. create new star particles
//    ===========================================================================================================
//    use OpenMP critical construct since both amr->Par->AddOneParticle() and amr->patch[0][lv][PID]->AddParticle()
//    will modify some global variables
//    --> note that the order of which thread calls amr->Par->AddOneParticle() is nondeterministic and may change from run to run
//        --> order of particles stored in the particle repository (i.e., their particle ID) may change from run to run
//        --> particle text file may change from run to run since it's dumped according to the order of particle ID
//    --> but it's not an issue since the actual data of each particle will not be affected
#     pragma omp critical
      {
//       4-1. add particles to the particle repository
         for (int p=0; p<NNewPar; p++)
            NewParID[p] = amr->Par->AddOneParticle( NewParAtt[p] );


//       4-2. add particles to the patch
         const real *PType = amr->Par->Type;
#        ifdef DEBUG_PARTICLE
//       do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//       may change after calling amr->Par->AddOneParticle()
         const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
         char Comment[100];
         sprintf( Comment, "%s", __FUNCTION__ );

         amr->patch[0][lv][PID]->AddParticle( NNewPar, NewParID, &amr->Par->NPar_Lv[lv],
                                              PType, ParPos, amr->Par->NPar_AcPlusInac, Comment );
#        else
         amr->patch[0][lv][PID]->AddParticle( NNewPar, NewParID, &amr->Par->NPar_Lv[lv], PType );
#        endif
      } // pragma omp critical

   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// free memory
   delete [] NewParAtt;
   delete [] NewParID;

   } // end of OpenMP parallel region


// get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

} // FUNCTION : SF_CreateStar_AGORA



#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )
