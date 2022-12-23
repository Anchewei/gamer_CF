#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Sink
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
                          const bool DetRandom, const bool UseMetal)
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
   if ( Idx_ParCreTime == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_ParCreTime is undefined !!\n" );
#  endif // #ifdef GAMER_DEBUG


// constant parameters
   const double dh             = amr->dh[lv];
   
   const int    AccCellNum     = 4;
   const double AccRadius      = AccCellNum*dh;
   const int    NGhost         = AccCellNum; // the number of ghost cell at each side
   const int    Size_Flu       = PS2 + 2*NGhost; // final cube size
   const int    Size_Flu_P1    = Size_Flu + 1; // for face-centered B field
   const int    Size_Pot       = Size_Flu; // for potential
   const int    NPG            = 1;

   const real   dv             = CUBE( dh );
   const int    FluSg          = amr->FluSg[lv];
   const real   Coeff_FreeFall = SQRT( (32.0*NEWTON_G)/(3.0*M_PI) );
// const real   GraConst       = ( OPT__GRA_P5_GRADIENT ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh);
   const real   GraConst       = ( false                ) ? -1.0/(12.0*dh) : -1.0/(2.0*dh); // P5 is NOT supported yet


// checking the value of accretion radius
   if (AccCellNum > PS1)
      Aux_Error( ERROR_INFO, "AccCellNum should be smaller than PATCH_SIZE !!" );

// start of OpenMP parallel region
#  pragma omp parallel
   {

// thread-private variables
#  ifdef OPENMP
   const int TID = omp_get_thread_num();
#  else
   const int TID = 0;
#  endif

   double x, y, z, VelX, VelY, VelZ;
   real   GasDens, GasDensFreeFall, StarMFrac, StarMass, GasMFracLeft;
   // real   (*fluid)[PS1][PS1][PS1]      = NULL; // fluid pointer, PS1: PATCH_SIZE

   const int MaxNewParPerPG = CUBE(PS2);
   real   (*NewParAtt)[PAR_NATT_TOTAL] = new real [MaxNewParPerPG][PAR_NATT_TOTAL];
   long    *NewParID                   = new long [MaxNewParPerPG];
   long    *NewParPID                  = new long [MaxNewParPerPG];

   int NNewPar;

// loop over all real patch groups
// use static schedule to ensure bitwise reproducibility when running with the same numbers of OpenMP threads and MPI ranks
// --> bitwise reproducibility will still break when running with different numbers of OpenMP threads and/or MPI ranks
//     unless both BITWISE_REPRODUCIBILITY and SF_CREATE_STAR_DET_RANDOM are enabled
#  pragma omp for schedule( static )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
//    skip non-leaf patches
      if ( amr->patch[0][lv][PID0]->son != -1 )  continue;

      real (*Flu_Array_F_In)[FLU_NIN][CUBE(Size_Flu)];
      real (*Mag_Array_F_In)[Size_Flu_P1*SQR(Size_Flu)];
      real (*Pot_Array_USG_F)[CUBE(Size_Pot)];
      real (*fluid)[FLU_NIN];
      real Corner_Array_F[3]; // the corner of the ghost zone

//    load the existing particles ID (the number)
      int   NParTot   = amr->Par->NPar_Active_AllRank;
      const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
      const real *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

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

      NNewPar = 0;
      for (int pk=NGhost; pk<PS2 + NGhost; pk++)
      for (int pj=NGhost; pj<PS2 + NGhost; pj++)
      for (int pi=NGhost; pi<PS2 + NGhost; pi++) // loop inside the patch group
      {
         x = Corner_Array_F[0] + pi*dh + dh*NGhost;
         y = Corner_Array_F[1] + pj*dh + dh*NGhost;
         z = Corner_Array_F[2] + pk*dh + dh*NGhost;

         const int t = IDX321( pi, pj, pk, Size_Flu, Size_Flu );
         for (int v=0; v<FLU_NIN; v++)    fluid[v] = Flu_Array_F_In[v][t];
         VelX = fluid[MOMX]/fluid[DENS];
         VelY = fluid[MOMY]/fluid[DENS];
         VelZ = fluid[MOMZ]/fluid[DENS];

//       1. Proximity Check + Density Threshold
//       ===========================================================================================================
         GasDens = fluid[DENS];
         if ( GasDens < GasDensThres )    continue;

         bool InsideAccRadius = false;
         bool NotPassDen      = false;

         for (int p=0; p<NParTot; p++)
         {
            real PCP[3]; // particle-cell relative position
            real D2Par; // particle-cell distance

            PCP[0] = x - ParPos[0][p];
            PCP[1] = y - ParPos[1][p];
            PCP[2] = z - ParPos[2][p];

            D2Par = SQRT(SQR(PCP[0])+SQR(PCP[1])+SQR(PCP[2]));
            if ( D2Par < AccRadius )
            {
               InsideAccRadius = true;
               break;
            }

            real PCV[3]; // particle-cell relative velocity
            real NPCP[3]; // normalized particle-cell relative position

            PCV[0] = VelX - ParVel[0][p];
            PCV[1] = VelY - ParVel[1][p];
            PCV[2] = VelZ - ParVel[2][p];

            NPCP[0] = PCP[0]/D2Par;
            NPCP[1] = PCP[1]/D2Par;
            NPCP[2] = PCP[2]/D2Par;

            GasDensFreeFall = SQR((1/Coeff_FreeFall)*(NPCP[0]*PCV[0] + NPCP[1]*PCV[1] + NPCP[2]*PCV[2])/D2Par); // Clarke et al. 2017, eqn (5)
            if ( GasDens < GasDensFreeFall )
            {
               NotPassDen = true;
               break;
            }
         } // NParTot

         if ( InsideAccRadius )               continue;
         if ( NotPassDen )                    continue;
         
//       2. Converging Flow Check
//       ===========================================================================================================
         real VelNeighbor[6]; // record the neighboring cell velocity [x+, x-, y+, y+, z+, z-]
         for (int NeighborID=0; NeighborID<6; NeighborID++)
         {  
            real dfluid[FLU_NIN]; // store the fluid in the adjacent cell
            if      (NeighborID == 0) int dt = IDX321(  1,  0,  0, Size_Flu, Size_Flu );
            else if (NeighborID == 1) int dt = IDX321( -1,  0,  0, Size_Flu, Size_Flu );
            else if (NeighborID == 2) int dt = IDX321(  0,  1,  0, Size_Flu, Size_Flu );
            else if (NeighborID == 3) int dt = IDX321(  0, -1,  0, Size_Flu, Size_Flu );
            else if (NeighborID == 4) int dt = IDX321(  0,  0,  1, Size_Flu, Size_Flu );
            else if (NeighborID == 5) int dt = IDX321(  0,  0, -1, Size_Flu, Size_Flu );

            int Neighbort = t + dt;
            for (int v=0; v<FLU_NIN; v++)    dfluid[v] = Flu_Array_F_In[v][Neighbort];

            if      ((NeighborID == 0) || (NeighborID == 1)) VelNeighbor[NeighborID] = dfluid[MOMX]/dfluid[DENS];
            else if ((NeighborID == 2) || (NeighborID == 3)) VelNeighbor[NeighborID] = dfluid[MOMY]/dfluid[DENS];
            else if ((NeighborID == 4) || (NeighborID == 5)) VelNeighbor[NeighborID] = dfluid[MOMZ]/dfluid[DENS];
         }

         real DivV = (VelNeighbor[0] + VelNeighbor[2] + VelNeighbor[4] - VelNeighbor[1] - VelNeighbor[3] - VelNeighbor[5]);
         if ( DivV > 0 )                       continue;

//       3. Gravitational Potential Minimum Check + Jeans Instability Check + Check for Bound State
//       ===========================================================================================================
         real Mtot = (real)0.0, MVel[3] = { (real)0.0, (real)0.0, (real)0.0}, MWvel[3]; // sum(mass_i), sum(mass_i*velocity_i), mass-weighted velocity
         for (int vk=0; vk<Size_Flu; vk++)
         for (int vj=0; vj<Size_Flu; vj++)
         for (int vi=0; vi<Size_Flu; vi++) // loop the all cells, to find the cells inside the control volumne (v)
         {  
            real (*vfluid)[FLU_NIN]; // store the fluid in the control volumne
            real vx = Corner_Array_F[0] + vi*dh;
            real vy = Corner_Array_F[1] + vj*dh;
            real vz = Corner_Array_F[2] + vk*dh;

            real D2CC = SQRT(SQR(vx - x)+SQR(vy - y)+SQR(vz - z)); // distance to the center cell
            if ( D2CC > AccRadius )                        continue; // check whether it is inside the control volume

            const int vt = IDX321( vi, vj, vk, Size_Flu, Size_Flu );
            for (int v=0; v<FLU_NIN; v++)    vfluid[v] = Flu_Array_F_In[v][vt];
            MVel[0] += vfluid[MOMX]*dv;
            MVel[1] += vfluid[MOMY]*dv;
            MVel[2] += vfluid[MOMZ]*dv;
            Mtot += vfluid[DENS]*dv;
         } // vi, vj, vk

         MWvel[0] = MVel[0]/Mtot;
         MWvel[1] = MVel[1]/Mtot;
         MWvel[2] = MVel[2]/Mtot; // COM velocity

         real CCEg = GasDens*Pot_Array_USG_F[t]; // Eg for the centered cell
         real Egtot = (real)0.0, Ethtot = (real)0.0, Emagtot = (real)0.0, Ekintot = (real)0.0;
         bool NotMiniEg      = false;
         for (int vk=0; vk<Size_Flu; vk++)
         for (int vj=0; vj<Size_Flu; vj++)
         for (int vi=0; vi<Size_Flu; vi++) // loop the all cells, to find the cells inside the control volumne (v)
         {
            real (*vfluid)[FLU_NIN]; // store the fluid in the control volumne
            real vx = Corner_Array_F[0] + vi*dh;
            real vy = Corner_Array_F[1] + vj*dh;
            real vz = Corner_Array_F[2] + vk*dh;

            real D2CC = SQRT(SQR(vx - x)+SQR(vy - y)+SQR(yz - z)); // distance to the center cell
            if ( D2CC > AccRadius )                        continue; // check whether it is inside the control volume

            const int vt = IDX321( vi, vj, vk, Size_Flu, Size_Flu );
            for (int v=0; v<FLU_NIN; v++)    vfluid[v] = Flu_Array_F_In[v][vt];

//          3.1 Gravitational Potential Minimum Check
            real vEg = vfluid[DENS]*Pot_Array_USG_F[vt]; // Eg for the current cell
            if ( vEg < CCEg )
            {
               NotMiniEg = true;
               break;
            }

//          3.2 Storing Egtot, Ethtot, Emagtot, Ekintot
            const bool CheckMinPres_No = false;
            real Pres, Cs2, vEmag=NULL_REAL;
            Egtot += vEg;

#           ifdef MHD
            vEmag = MHD_GetCellCenteredBEnergy( Mag_Array_F_In[MAGX],
                                                Mag_Array_F_In[MAGY],
                                                Mag_Array_F_In[MAGZ],
                                                Size_Flu, Size_Flu, Size_Flu, vi, vj, vk );
            Emagtot += vEmag;
#           endif

            Pres = Hydro_Con2Pres( vfluid[DENS], vfluid[MOMX], vfluid[MOMY], vfluid[MOMZ], vfluid[ENGY],
                                   vfluid+NCOMP_FLUID, CheckMinPres_No, NULL_REAL, vEmag,
                                   EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
            Cs2  = EoS_DensPres2CSqr_CPUPtr( vfluid[DENS], Pres, vfluid+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
            Ethtot += 0.5*vfluid[DENS]*Cs2;

            Ekintot += 0.5*vfluid[DENS]*( SQR(vfluid[MOMX]/vfluid[DENS] - MWvel[0]) + SQR(vfluid[MOMY]/vfluid[DENS] - MWvel[1]) + SQR(vfluid[MOMZ]/vfluid[DENS] - MWvel[2]));
         } // vi, vj, vk

         if ( NotMiniEg )                                   continue;
         if ( FABS(Egtot) < 2*Ethtot)                       continue;
         if (( Egtot + Ethtot + Ekintot + Emagtot ) > 0)    continue;

//       4. store the information of new star particles
//       --> we will not create these new particles until looping over all cells in a patch in order to reduce
//           the OpenMP synchronization overhead
//       ===========================================================================================================
#        ifdef GAMER_DEBUG
         if ( NNewPar >= MaxNewParPerPG )
            Aux_Error( ERROR_INFO, "NNewPar (%d) >= MaxNewParPerPG (%d) !!\n", NNewPar, MaxNewParPerPG );
#        endif

         NewParAtt[NNewPar][PAR_MASS] = (GasDens - GasDensThres)*dv;
         NewParAtt[NNewPar][PAR_POSX] = x;
         NewParAtt[NNewPar][PAR_POSY] = y;
         NewParAtt[NNewPar][PAR_POSZ] = z;
         NewParAtt[NNewPar][PAR_VELX] = VelX;
         NewParAtt[NNewPar][PAR_VELY] = VelY;
         NewParAtt[NNewPar][PAR_VELZ] = VelZ;
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
            real PotNeighbor[6]; // record the neighboring cell potential [x+, x-, y+, y+, z+, z-]
            for (int NeighborID=0; NeighborID<6; NeighborID++)
            {  
               if      (NeighborID == 0) int dt = IDX321(  1,  0,  0, Size_Flu, Size_Flu );
               else if (NeighborID == 1) int dt = IDX321( -1,  0,  0, Size_Flu, Size_Flu );
               else if (NeighborID == 2) int dt = IDX321(  0,  1,  0, Size_Flu, Size_Flu );
               else if (NeighborID == 3) int dt = IDX321(  0, -1,  0, Size_Flu, Size_Flu );
               else if (NeighborID == 4) int dt = IDX321(  0,  0,  1, Size_Flu, Size_Flu );
               else if (NeighborID == 5) int dt = IDX321(  0,  0, -1, Size_Flu, Size_Flu );

               int Neighbort = t + dt;
               PotNeighbor[NeighborID] = Pot_Array_USG_F[Neighbort];
            }
            GasAcc[0] += PotNeighbor[0] - PotNeighbor[1];
            GasAcc[1] += PotNeighbor[2] - PotNeighbor[3];
            GasAcc[2] += PotNeighbor[4] - PotNeighbor[5];
         }

         NewParAtt[NNewPar][PAR_ACCX] = GasAcc[0];
         NewParAtt[NNewPar][PAR_ACCY] = GasAcc[1];
         NewParAtt[NNewPar][PAR_ACCZ] = GasAcc[2];
#        endif // ifdef STORE_PAR_ACC

         NewParAtt[NNewPar][Idx_ParCreTime  ] = TimeNew;

//       5. remove the gas that has been converted to stars
//       ===========================================================================================================
         GasMFracLeft = (real) 1.0 - (GasDensThres/GasDens);
         int PGi = pi - NGhost;
         int PGj = pj - NGhost;
         int PGk = pk - NGhost; // the cell id inside patch group
         
         // determine the current patch where the sink is
         int LocalID = 0;
         if      ((PGi < PS1) && (PGj < PS1) && (PGk < PS1))     LocalID = 0;
         else if ((PGi >= PS1) && (PGj < PS1) && (PGk < PS1))    LocalID = 1;
         else if ((PGi < PS1) && (PGj >= PS1) && (PGk < PS1))    LocalID = 2;
         else if ((PGi < PS1) && (PGj < PS1) && (PGk >= PS1))    LocalID = 3;
         else if ((PGi >= PS1) && (PGj >= PS1) && (PGk < PS1))   LocalID = 4;
         else if ((PGi < PS1) && (PGj >= PS1) && (PGk >= PS1))   LocalID = 5;
         else if ((PGi >= PS1) && (PGj < PS1) && (PGk >= PS1))   LocalID = 6;
         else if ((PGi >= PS1) && (PGj >= PS1) && (PGk >= PS1))  LocalID = 7;

         const int Disp_i = TABLE_02( LocalID, 'x', 0, PS1 );
         const int Disp_j = TABLE_02( LocalID, 'y', 0, PS1 );
         const int Disp_k = TABLE_02( LocalID, 'z', 0, PS1 );
         const int PID = PID0 + LocalID;
         NewParPID[NNewPar] = PID;

         for (int v=0; v<NCOMP_TOTAL; v++)
         amr->patch[FluSg][lv][PID]->fluid[v][PGk - Disp_k][PGj - Disp_j][PGi - Disp_i] *= GasMFracLeft;

         NNewPar ++;
      } // pi, pj, pk

//    6. create new star particles
//    ===========================================================================================================
//    use OpenMP critical construct since both amr->Par->AddOneParticle() and amr->patch[0][lv][PID]->AddParticle()
//    will modify some global variables
//    --> note that the order of which thread calls amr->Par->AddOneParticle() is nondeterministic and may change from run to run
//        --> order of particles stored in the particle repository (i.e., their particle ID) may change from run to run
//        --> particle text file may change from run to run since it's dumped according to the order of particle ID
//    --> but it's not an issue since the actual data of each particle will not be affected
#     pragma omp critical
      {
//       6-1. add particles to the particle repository
         for (int p=0; p<NNewPar; p++)
            NewParID[p] = amr->Par->AddOneParticle( NewParAtt[p] );


//       6-2. add particles to the patch
         const real *PType = amr->Par->Type;
#        ifdef DEBUG_PARTICLE
//       do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//       may change after calling amr->Par->AddOneParticle()
         const real *NewParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
         char Comment[100];
         sprintf( Comment, "%s", __FUNCTION__ );
         
         for (int p=0; p<NNewPar; p++) // since the particles can be in different PID, we add them one by one
         amr->patch[0][lv][NewParPID[p]]->AddParticle( 1, *NewParID[p], &amr->Par->NPar_Lv[lv],
                                                       PType, NewParPos, amr->Par->NPar_AcPlusInac, Comment );

#        else
         for (int p=0; p<NNewPar; p++) // since the particles can be in different PID, we add them one by one
         amr->patch[0][lv][NewParPID[p]]->AddParticle( 1, *NewParID[p], &amr->Par->NPar_Lv[lv], PType );
#        endif
      } // pragma omp critical

   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


// free memory
   delete [] NewParAtt;
   delete [] NewParID;

   } // end of OpenMP parallel region

// get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

} // FUNCTION : SF_CreateStar_Sink

#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )