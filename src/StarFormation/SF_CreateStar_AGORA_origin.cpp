#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  FindLocalPID
// Description :  Find the local PID relative to PID0 wihtin a patch group + ghost zone (size given by NGhost)
// Note        :  1. Here, we assume each patch consists of 8 patches.
//                2. The full size of a patch group + ghost zone = CUBE(PS2+2*NGhost)
// Parameter   :  i, j, k           : The cell index in a patch group + ghost zone
//                PS1               : The size of one patch
//                NGosht            : The size of ghost zone padded around the patch group 
// Return      :  The local PID relative to PID0
//-------------------------------------------------------------------------------------------------------

int FindLocalPID(int pi, int pj, int pk, int &PGi, int &PGj, int &PGk, int NGhost)
{
   PGi = pi - NGhost;
   PGj = pj - NGhost;
   PGk = pk - NGhost; // the cell id inside patch group
   
   // determine the current patch where the current cell is
   int CellPID;
   if      ((PGi < PS1) && (PGj < PS1) && (PGk < PS1))     CellPID = 0;
   else if ((PGi >= PS1) && (PGj < PS1) && (PGk < PS1))    CellPID = 1;
   else if ((PGi < PS1) && (PGj >= PS1) && (PGk < PS1))    CellPID = 2;
   else if ((PGi < PS1) && (PGj < PS1) && (PGk >= PS1))    CellPID = 3;
   else if ((PGi >= PS1) && (PGj >= PS1) && (PGk < PS1))   CellPID = 4;
   else if ((PGi < PS1) && (PGj >= PS1) && (PGk >= PS1))   CellPID = 5;
   else if ((PGi >= PS1) && (PGj < PS1) && (PGk >= PS1))   CellPID = 6;
   else if ((PGi >= PS1) && (PGj >= PS1) && (PGk >= PS1))  CellPID = 7;
   return CellPID;
}

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

#  ifdef MY_DEBUG
   const char  FileName[] = "Record__Par_debug";
   FILE *File = fopen( FileName, "a" );
#  endif

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

   real x, y, z, VelX, VelY, VelZ, vx, vy, vz;
   real GasDens, GasDensFreeFall, GasMFracLeft;
   real fluid[FLU_NIN], dfluid[FLU_NIN], vfluid[FLU_NIN];
   real Corner_Array_F[3]; // the corner of the ghost zone
   real PCP[3], D2Par, PCV[3]; // particle-cell relative position, distance, relative velocity
   real D2CC; // distance to the center cell
   real NPCP[3]; // normalized particle-cell relative position
   real VelNeighbor[6]; // record the neighboring cell velocity [x+, x-, y+, y+, z+, z-]
   real vEg, Pres, Cs2, vEmag=NULL_REAL;
   real PotNeighbor[6]; // record the neighboring cell potential [x+, x-, y+, y+, z+, z-]


   const int MaxNewParPerPG = CUBE(PS2);
   real   (*NewParAtt)[PAR_NATT_TOTAL]        = new real [MaxNewParPerPG][PAR_NATT_TOTAL];
   long    *NewParID                          = new long [MaxNewParPerPG];
   long    *NewParPID                         = new long [MaxNewParPerPG];
   real   (*Flu_Array_F_In)[CUBE(Size_Flu)]   = new real [FLU_NIN][CUBE(Size_Flu)];
   real   (*Mag_Array_F_In)                   = new real [Size_Flu_P1*SQR(Size_Flu)];
   real   (*Pot_Array_USG_F)                  = new real [CUBE(Size_Pot)];

   int NNewPar, LocalID, delta_t, PGi, PGj, PGk;

//    load the existing particles ID (their number)
   // const int   NParTot   = amr->Par->NPar_Active_AllRank;
   // const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   // const real *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };

   const bool TimingSendPar_Yes = true;
   const bool JustCountNPar_No  = false;
   const bool PredictPos_No     = false;
   const bool SibBufPatch_Yes   = true;
   const bool FaSibBufPatch_No  = false;
   const long ParAttBitIdx_In   = _PAR_TOTAL;
   Par_CollectParticle2OneLevel( lv, ParAttBitIdx_In, PredictPos_No, TimeNew, SibBufPatch_Yes, FaSibBufPatch_No,
                                 JustCountNPar_No, TimingSendPar_Yes );

// get the sibling index differences along different directions
   int NSibPID_Delta[26], *SibPID_Delta[26];

   TABLE_GetSibPID_Delta( NSibPID_Delta, SibPID_Delta );

// loop over all real patch groups
// use static schedule to ensure bitwise reproducibility when running with the same numbers of OpenMP threads and MPI ranks
// --> bitwise reproducibility will still break when running with different numbers of OpenMP threads and/or MPI ranks
//     unless both BITWISE_REPRODUCIBILITY and SF_CREATE_STAR_DET_RANDOM are enabled
#  pragma omp for schedule( static )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
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
      real *Mag_Array = Mag_Array_F_In[0];
#     else
      real *Mag_Array = NULL;
#     endif
      Prepare_PatchData( lv, TimeNew, Flu_Array_F_In[0], Mag_Array,
                        NGhost, NPG, &PID0, _TOTAL, _MAG,
                        OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                        OPT__BC_FLU, BC_POT_NONE, MinDens,    MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency );

#     ifdef UNSPLIT_GRAVITY
//    prepare the potential array
      if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
      Prepare_PatchData( lv, TimeNew, Pot_Array_USG_F, NULL,
                        NGhost, NPG, &PID0, _POTE, _NONE,
                        OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                        OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

//    prepare the corner array
      for (int d=0; d<3; d++)    Corner_Array_F[d] = amr->patch[0][lv][PID0]->EdgeL[d] + 0.5*dh - dh*NGhost;

#     endif // #ifdef UNSPLIT_GRAVITY


//    prepare the sibling patches id for this patch group
//    collect patch information
      const int NNearbyPatchMax = 64;  // maximum number of neaby patches of a patch group (including 8 local patches)
      int Nearby_PID_List[NNearbyPatchMax], NNearbyPatch, SibPID0_List[26];

//    get nearby patches
      NNearbyPatch = 0;

//    local patches
      for (int PID=PID0; PID<PID0+8; PID++)  Nearby_PID_List[ NNearbyPatch ++ ] = PID;

//    sibling patches
      TABLE_GetSibPID_Based( lv, PID0, SibPID0_List );

//###OPTIMIZATION: skip sibling patches if the maximum feedback radius is zero
      for (int s=0; s<26; s++)
      {
         const int SibPID0 = SibPID0_List[s];   // first target patch in the sibling patch group

//       only consider leaf patches on FB_LEVEL (including both real and buffer patches)
         if ( SibPID0 >= 0 )
         for (int c=0; c<NSibPID_Delta[s]; c++)
         {
            const int SibPID = SibPID0 + SibPID_Delta[s][c];
            Nearby_PID_List[ NNearbyPatch ++ ] = SibPID;
         }
      }


      NNewPar = 0;
      for (int pk=NGhost; pk<PS2 + NGhost; pk++)
      for (int pj=NGhost; pj<PS2 + NGhost; pj++)
      for (int pi=NGhost; pi<PS2 + NGhost; pi++) // loop inside the patch group
      {  
         LocalID = FindLocalPID(pi, pj, pk, PGi, PGj, PGk, NGhost);

         const int Disp_i = TABLE_02( LocalID, 'x', 0, PS1 ); // the cell index within PID
         const int Disp_j = TABLE_02( LocalID, 'y', 0, PS1 );
         const int Disp_k = TABLE_02( LocalID, 'z', 0, PS1 );
         const int PID = PID0 + LocalID; // record the current PID

//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         x = Corner_Array_F[0] + pi*dh;
         y = Corner_Array_F[1] + pj*dh;
         z = Corner_Array_F[2] + pk*dh;

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

         real  *ParAtt_Local[PAR_NATT_TOTAL];
         int    NParMax   = -1;

         for (int v=0; v<PAR_NATT_TOTAL; v++)   ParAtt_Local[v] = NULL;

         for (int t=0; t<NNearbyPatch; t++)
         {
            const int PPID = Nearby_PID_List[t]; // check the particle number for this patch (PPID)

         // check both NPar and NPar_Copy (NPar_Copy may be -1, which is fine)
            NParMax = MAX( NParMax, amr->patch[0][lv][PPID]->NPar      );
            NParMax = MAX( NParMax, amr->patch[0][lv][PPID]->NPar_Copy );
         }

         if ( NParMax > 0 )
         {
            for (int v=0; v<PAR_NATT_TOTAL; v++)
               if ( ParAttBitIdx_In & BIDX(v) )    ParAtt_Local[v] = new real [NParMax];
         }

         // iterate over all nearby patches of the target patch group
         for (int t=0; t<NNearbyPatch; t++)
         {
            const int PPID = Nearby_PID_List[t];
            long  *ParList = NULL;
            int    NPar;
            bool   UseParAttCopy;

//          the check "son == -1" is actually useless for now since we fix LEVEL == MAX_LEVEL
            if ( amr->patch[0][lv][PPID]->son == -1  &&  PPID < amr->NPatchComma[lv][1] )
            {
               NPar          = amr->patch[0][lv][PPID]->NPar;
               ParList       = amr->patch[0][lv][PPID]->ParList;
               UseParAttCopy = false;

#              ifdef DEBUG_PARTICLE
               if ( amr->patch[0][lv][PPID]->NPar_Copy != -1 )
                  Aux_Error( ERROR_INFO, "lv %d, PPID %d, NPar_Copy = %d != -1 !!\n",
                           lv, PPID, amr->patch[0][lv][PPID]->NPar_Copy );
#              endif
            }

            else
            {
//             note that amr->patch[0][lv][PPID]->NPar>0 is still possible
               NPar          = amr->patch[0][lv][PPID]->NPar_Copy;
#              ifdef LOAD_BALANCE
               ParList       = NULL;
               UseParAttCopy = true;
#              else
               ParList       = amr->patch[0][lv][PPID]->ParList_Copy;
               UseParAttCopy = false;
#              endif
            } // if ( amr->patch[0][lv][PPID]->son == -1  &&  PPID < amr->NPatchComma[lv][1] ) ... else ...

#           ifdef LOAD_BALANCE
            if ( UseParAttCopy ) {
               for (int v=0; v<PAR_NATT_TOTAL; v++) {
                  if ( ParAttBitIdx_In & BIDX(v) ) {

#                 ifdef DEBUG_PARTICLE
                  if ( NPar > 0  &&  amr->patch[0][lv][PPID]->ParAtt_Copy[v] == NULL )
                     Aux_Error( ERROR_INFO, "ParAtt_Copy == NULL for NPar (%d) > 0 (lv %d, PPID %d, v %d) !!\n",
                                 NPar, lv, PPID, v );
#                 endif

                  for (int p=0; p<NPar; p++)
                     ParAtt_Local[v][p] = amr->patch[0][lv][PPID]->ParAtt_Copy[v][p];
            }}}

            else
#           endif // #ifdef LOAD_BALANCE
            {
#              ifdef DEBUG_PARTICLE
               if ( NPar > 0  &&  ParList == NULL )
                  Aux_Error( ERROR_INFO, "ParList == NULL for NPar (%d) > 0 (lv %d, PPID %d) !!\n",
                              NPar, lv, PPID );
#              endif

               for (int v=0; v<PAR_NATT_TOTAL; v++) {
                  if ( ParAttBitIdx_In & BIDX(v) )
                     for (int p=0; p<NPar; p++)
                        ParAtt_Local[v][p] = amr->Par->Attribute[v][ ParList[p] ];
               }
            } // if ( UseParAttCopy ) ... else ...

            for (int p=0; p<NPar; p++) // loop over all nearby particles
            {
               PCP[0] = x - ParAtt_Local[PAR_POSX][p];
               PCP[1] = y - ParAtt_Local[PAR_POSY][p];
               PCP[2] = z - ParAtt_Local[PAR_POSZ][p];
               D2Par = SQRT(SQR(PCP[0])+SQR(PCP[1])+SQR(PCP[2]));
               if ( D2Par < AccRadius )
               {
                  InsideAccRadius = true;
                  break;
               }

               PCV[0] = VelX - ParAtt_Local[PAR_VELX][p];
               PCV[1] = VelY - ParAtt_Local[PAR_VELY][p];
               PCV[2] = VelZ - ParAtt_Local[PAR_VELZ][p];

               NPCP[0] = PCP[0]/D2Par;
               NPCP[1] = PCP[1]/D2Par;
               NPCP[2] = PCP[2]/D2Par;

               GasDensFreeFall = SQR((1/Coeff_FreeFall)*(NPCP[0]*PCV[0] + NPCP[1]*PCV[1] + NPCP[2]*PCV[2])/D2Par); // Clarke et al. 2017, eqn (5)
               if ( GasDens < GasDensFreeFall )
               {
                  NotPassDen = true;
                  break;
               }

            } // for (int p=0; p<NPar; p++) 

            if ( InsideAccRadius )           break;
            if ( NotPassDen )                break;
         } // for (int t=0; t<NNearbyPatch; t++)

         if ( InsideAccRadius )               continue;
         if ( NotPassDen )                    continue;

         for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParAtt_Local[v];
         
//       2. Converging Flow Check
//       ===========================================================================================================
         for (int NeighborID=0; NeighborID<6; NeighborID++)
         {  
            if      (NeighborID == 0) delta_t = IDX321(  1,  0,  0, Size_Flu, Size_Flu );
            else if (NeighborID == 1) delta_t = IDX321( -1,  0,  0, Size_Flu, Size_Flu );
            else if (NeighborID == 2) delta_t = IDX321(  0,  1,  0, Size_Flu, Size_Flu );
            else if (NeighborID == 3) delta_t = IDX321(  0, -1,  0, Size_Flu, Size_Flu );
            else if (NeighborID == 4) delta_t = IDX321(  0,  0,  1, Size_Flu, Size_Flu );
            else if (NeighborID == 5) delta_t = IDX321(  0,  0, -1, Size_Flu, Size_Flu );

            const int Neighbort = t + delta_t;
            for (int v=0; v<FLU_NIN; v++)    dfluid[v] = Flu_Array_F_In[v][Neighbort];

            if      ((NeighborID == 0) || (NeighborID == 1)) VelNeighbor[NeighborID] = dfluid[MOMX]/dfluid[DENS];
            else if ((NeighborID == 2) || (NeighborID == 3)) VelNeighbor[NeighborID] = dfluid[MOMY]/dfluid[DENS];
            else if ((NeighborID == 4) || (NeighborID == 5)) VelNeighbor[NeighborID] = dfluid[MOMZ]/dfluid[DENS];
         } // for (int NeighborID=0; NeighborID<6; NeighborID++)

         if ( (VelNeighbor[0] - VelNeighbor[1]) >= 0 )                       continue;
         if ( (VelNeighbor[2] - VelNeighbor[3]) >= 0 )                       continue;
         if ( (VelNeighbor[4] - VelNeighbor[5]) >= 0 )                       continue;

//       3. Gravitational Potential Minimum Check + Jeans Instability Check + Check for Bound State
//       ===========================================================================================================
         real Mtot = (real)0.0, MVel[3] = { (real)0.0, (real)0.0, (real)0.0}, MWvel[3]; // sum(mass_i), sum(mass_i*velocity_i), mass-weighted velocity
         real phi000 = (real)0.0;
// #        ifdef MY_DEBUG
//          int count = 0;
// #        endif
         for (int vk=pk-AccCellNum; vk<=pk+AccCellNum; vk++)
         for (int vj=pj-AccCellNum; vj<=pj+AccCellNum; vj++)
         for (int vi=pi-AccCellNum; vi<=pi+AccCellNum; vi++) // loop the nearby cells, to find the cells inside the control volumne (v)
         {
            vx = Corner_Array_F[0] + vi*dh;
            vy = Corner_Array_F[1] + vj*dh;
            vz = Corner_Array_F[2] + vk*dh;

            D2CC = SQRT(SQR(vx - x)+SQR(vy - y)+SQR(vz - z)); // distance to the center cell
            if ( D2CC > AccRadius )                 continue; // check whether it is inside the control volume

            const int vt = IDX321( vi, vj, vk, Size_Flu, Size_Flu );
            for (int v=0; v<FLU_NIN; v++)    vfluid[v] = Flu_Array_F_In[v][vt];
            MVel[0] += vfluid[MOMX]*dv;
            MVel[1] += vfluid[MOMY]*dv;
            MVel[2] += vfluid[MOMZ]*dv;
            Mtot += vfluid[DENS]*dv;

            if ( D2CC != 0.0 )        phi000 += -NEWTON_G*vfluid[DENS]*dv/D2CC; // potential
// #           ifdef MY_DEBUG
//             if ((PID == 31) and (pi==19) and (pj==19)and ( pk==19))
//             {
//             count++;
//             fprintf( File, "%7.4e %d %d%d%d %d %7.4e %7.4e %7.4e",  TimeNew, PID, vi, vj, vk, count,
//                      phi000, Mtot, MVel[0]);
//             fprintf( File, "\n" );
//             }
//             // fprintf( File, "'%13.7e %13.7e %13.7e',",  x, y, z);
//             // fprintf( File, "\n" );
// #           endif
         } // vi, vj, vk

         MWvel[0] = MVel[0]/Mtot;
         MWvel[1] = MVel[1]/Mtot;
         MWvel[2] = MVel[2]/Mtot; // COM velocity

         // Calculate for the surrounding cells and totoal energy
         real Egtot = (real)0.0;
         real vifluid[FLU_NIN], vjfluid[FLU_NIN];
         bool NotMiniPot      = false; // do the check during the process
         for (int vki=pk-AccCellNum; vki<=pk+AccCellNum; vki++)
         for (int vji=pj-AccCellNum; vji<=pj+AccCellNum; vji++)
         for (int vii=pi-AccCellNum; vii<=pi+AccCellNum; vii++) // loop the nearby cells, to find the cells inside the control volumne (v)
         {
            real vxi, vyi, vzi, vxj, vyj, vzj;
            real phiijk = (real)0.0;
            vxi = Corner_Array_F[0] + vii*dh;
            vyi = Corner_Array_F[1] + vji*dh;
            vzi = Corner_Array_F[2] + vki*dh;

            real D2CCi = SQRT(SQR(vxi - x)+SQR(vyi - y)+SQR(vzi - z)); // distance to the center cell
            if ( D2CCi > AccRadius )                 continue; // check whether it is inside the control volume
            const int vti = IDX321( vii, vji, vki, Size_Flu, Size_Flu );
            for (int v=0; v<FLU_NIN; v++)    vifluid[v] = Flu_Array_F_In[v][vti];

            for (int vkj=pk-AccCellNum; vkj<=pk+AccCellNum; vkj++)
            for (int vjj=pj-AccCellNum; vjj<=pj+AccCellNum; vjj++)
            for (int vij=pi-AccCellNum; vij<=pi+AccCellNum; vij++) // loop the nearby cells, to find the cells inside the control volumne (v)
            {
               vxj = Corner_Array_F[0] + vij*dh;
               vyj = Corner_Array_F[1] + vjj*dh;
               vzj = Corner_Array_F[2] + vkj*dh;

               real rij = SQRT(SQR(vxi - vxj)+SQR(vyi - vyj)+SQR(vzi - vzj));
               if ( rij == 0.0 )                        continue;

               real D2CCj = SQRT(SQR(vxj - x)+SQR(vyj - y)+SQR(vzj - z)); // distance to the center cell
               if ( D2CCj > AccRadius )                 continue; // check whether it is inside the control volume
               const int vtj = IDX321( vij, vjj, vkj, Size_Flu, Size_Flu );
               for (int v=0; v<FLU_NIN; v++)    vjfluid[v] = Flu_Array_F_In[v][vtj];

               phiijk += -NEWTON_G*vjfluid[DENS]*dv/rij; // potential
            } // vij, vjj, vkj

//          3.1 Gravitational Potential Minimum Check
            if ( phi000 != MIN( phi000, phiijk ) )
            {
               NotMiniPot = true;
               break;
            }

            Egtot += 0.5*vifluid[DENS]*dv*phiijk;
         } // vii, vji, vki
// #        ifdef MY_DEBUG
//          const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
//          real dx, dy, dz;
//          dx = x - BoxCenter[0];
//          dy = y - BoxCenter[1];
//          dz = z - BoxCenter[2];
//          real DOC = SQRT(SQR(dx)+SQR(dy)+SQR(dz));
//          if (DOC <= 0.5*AccRadius)
//          {
//          fprintf( File, "%7.4e %d %d%d%d %d",  TimeNew, PID, pi, pj, pk,
//                   NotMiniPot);
//          fprintf( File, "\n" );
//          }
//          // fprintf( File, "'%13.7e %13.7e %13.7e',",  x, y, z);
//          // fprintf( File, "\n" );
// #        endif
         if ( NotMiniPot )                                   continue;

         real Ethtot = (real)0.0, Emagtot = (real)0.0, Ekintot = (real)0.0;
         for (int vk=pk-AccCellNum; vk<=pk+AccCellNum; vk++)
         for (int vj=pj-AccCellNum; vj<=pj+AccCellNum; vj++)
         for (int vi=pi-AccCellNum; vi<=pi+AccCellNum; vi++) // loop the nearby cells, to find the cells inside the control volumne (v)
         {
            vx = Corner_Array_F[0] + vi*dh;
            vy = Corner_Array_F[1] + vj*dh;
            vz = Corner_Array_F[2] + vk*dh;

            D2CC = SQRT(SQR(vx - x)+SQR(vy - y)+SQR(vz - z)); // distance to the center cell
            if ( D2CC > AccRadius )                 continue; // check whether it is inside the control volume

            const int vt = IDX321( vi, vj, vk, Size_Flu, Size_Flu );
            for (int v=0; v<FLU_NIN; v++)    vfluid[v] = Flu_Array_F_In[v][vt];

//          3.2 Storing Egtot, Ethtot, Emagtot, Ekintot
            const bool CheckMinPres_No = false;

#           ifdef MHD
            vEmag = MHD_GetCellCenteredBEnergy( Mag_Array_F_In[MAGX],
                                                Mag_Array_F_In[MAGY],
                                                Mag_Array_F_In[MAGZ],
                                                Size_Flu, Size_Flu, Size_Flu, vi, vj, vk );
            Emagtot += vEmag*dv;
#           endif

            Pres = Hydro_Con2Pres( vfluid[DENS], vfluid[MOMX], vfluid[MOMY], vfluid[MOMZ], vfluid[ENGY],
                                   vfluid+NCOMP_FLUID, CheckMinPres_No, NULL_REAL, vEmag,
                                   EoS_DensEint2Pres_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
            Cs2  = EoS_DensPres2CSqr_CPUPtr( vfluid[DENS], Pres, vfluid+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
            Ethtot += 0.5*vfluid[DENS]*dv*Cs2;

            Ekintot += 0.5*vfluid[DENS]*dv*( SQR(vfluid[MOMX]/vfluid[DENS] - MWvel[0]) + SQR(vfluid[MOMY]/vfluid[DENS] - MWvel[1]) + SQR(vfluid[MOMZ]/vfluid[DENS] - MWvel[2]));
         } // vi, vj, vk

         if ( FABS(Egtot) <= 2*Ethtot)                       continue;
         if (( Egtot + Ethtot + Ekintot + Emagtot ) >= 0)    continue;
#        ifdef MY_DEBUG
         fprintf( File, "'%d %13.7e %13.7e %13.7e %13.7e',",  PID, x, y, z, phi000);
         fprintf( File, "\n" );
#        endif
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
            for (int NeighborID=0; NeighborID<6; NeighborID++)
            {  
               if      (NeighborID == 0) delta_t = IDX321(  1,  0,  0, Size_Flu, Size_Flu );
               else if (NeighborID == 1) delta_t = IDX321( -1,  0,  0, Size_Flu, Size_Flu );
               else if (NeighborID == 2) delta_t = IDX321(  0,  1,  0, Size_Flu, Size_Flu );
               else if (NeighborID == 3) delta_t = IDX321(  0, -1,  0, Size_Flu, Size_Flu );
               else if (NeighborID == 4) delta_t = IDX321(  0,  0,  1, Size_Flu, Size_Flu );
               else if (NeighborID == 5) delta_t = IDX321(  0,  0, -1, Size_Flu, Size_Flu );

               const int Neighbort = t + delta_t;
               PotNeighbor[NeighborID] = Pot_Array_USG_F[Neighbort];
            }
            GasAcc[0] += GraConst*(PotNeighbor[0] - PotNeighbor[1]);
            GasAcc[1] += GraConst*(PotNeighbor[2] - PotNeighbor[3]);
            GasAcc[2] += GraConst*(PotNeighbor[4] - PotNeighbor[5]);
         }

         NewParAtt[NNewPar][PAR_ACCX] = GasAcc[0];
         NewParAtt[NNewPar][PAR_ACCY] = GasAcc[1];
         NewParAtt[NNewPar][PAR_ACCZ] = GasAcc[2];
#        endif // ifdef STORE_PAR_ACC

         NewParAtt[NNewPar][Idx_ParCreTime  ] = TimeNew;

//       5. remove the gas that has been converted to stars
//       ===========================================================================================================
         GasMFracLeft = (real) 1.0 - (GasDensThres/GasDens);
         NewParPID[NNewPar] = PID;

         for (int v=0; v<NCOMP_TOTAL; v++)
         amr->patch[FluSg][lv][PID]->fluid[v][PGk - Disp_k][PGj - Disp_j][PGi - Disp_i] *= GasMFracLeft;
         
         NNewPar ++;
      } // pi, pj, pk

   // 6. create new star particles
   // ===========================================================================================================
   // use OpenMP critical construct since both amr->Par->AddOneParticle() and amr->patch[0][lv][PID]->AddParticle()
   // will modify some global variables
   // --> note that the order of which thread calls amr->Par->AddOneParticle() is nondeterministic and may change from run to run
   //     --> order of particles stored in the particle repository (i.e., their particle ID) may change from run to run
   //     --> particle text file may change from run to run since it's dumped according to the order of particle ID
   // --> but it's not an issue since the actual data of each particle will not be affected
#     pragma omp critical
      {
//       6-1. add particles to the particle repository
         for (int p=0; p<NNewPar; p++)
            NewParID[p] = amr->Par->AddOneParticle( NewParAtt[p] );


//       6-2. add particles to the patch
         const real *PType = amr->Par->Type;
         int ParInPatch;
         for (int PID=PID0; PID<PID0+8; PID++)
         {
            long    *ParIDInPatch      = new long [MaxNewParPerPG]; // ParID in the current patch
            ParInPatch = 0;
            for (int p=0; p<NNewPar; p++)
            {
               if (NewParPID[p] == PID) 
               {
                  ParIDInPatch[ParInPatch] = NewParID[p];
                  ParInPatch ++;
               } // if (NewParPID[p] == PID) 
            } // for (int p=0; p<NNewPar; p++)

            if ( ParInPatch == 0 )                        continue;

#           ifdef DEBUG_PARTICLE
//          do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//          may change after calling amr->Par->AddOneParticle()
            const real *NewParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
            char Comment[100];
            sprintf( Comment, "%s", __FUNCTION__ );
            
            amr->patch[0][lv][PID]->AddParticle(ParInPatch, ParIDInPatch, &amr->Par->NPar_Lv[lv],
                                                         PType, NewParPos, amr->Par->NPar_AcPlusInac, Comment );

#           else
            amr->patch[0][lv][PID]->AddParticle( ParInPatch, ParIDInPatch, &amr->Par->NPar_Lv[lv], PType );
#           endif

            delete [] ParIDInPatch;

         } // for (int PID=PID0; PID<PID0+8; PID++)
      } // pragma omp critical
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

#  ifdef MY_DEBUG
   fprintf( File, "Step finished");
   fprintf( File, "\n" );
   fclose( File );
#  endif


// free memory
   delete [] NewParAtt;
   delete [] NewParID;
   delete [] NewParPID;
   delete [] Flu_Array_F_In;
   delete [] Mag_Array_F_In;
   delete [] Pot_Array_USG_F;
   Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch_Yes, FaSibBufPatch_No );
   } // end of OpenMP parallel region

// get the total number of active particles in all MPI ranks
   MPI_Allreduce( &amr->Par->NPar_Active, &amr->Par->NPar_Active_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

} // FUNCTION : SF_CreateStar_Sink

#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )