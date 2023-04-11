#include "GAMER.h"

#ifdef FEEDBACK

extern double AccGasDensThres;

// function pointers to be set by FB_Init_Plummer()
extern int (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                           const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                           real (*Fluid)[FB_NXT][FB_NXT][FB_NXT], const double EdgeL[], const double dh, bool CoarseFine[],
                           const int TID, RandomNumber_t *RNG );
extern void (*FB_End_User_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Plummer
// Description :  Example of explosion and mass accretion feedbacks
//                --> Just an example and is by no means to be physically correct
//
// Note        :  1. Input and output fluid and particle data are stored in Fluid[] and ParAtt[], respectively
//                   --> This function is responsible for updating gas and particles within
//                       ** FB_GHOST_SIZE <= cell indices i,j,k < FB_GHOST_SIZE+PS2 **
//                   --> Updating gas and particles outside this range is fine but will have no effect at all
//                2. Must use ParSortID[] to access ParAtt[]
//                   --> ParAtt[PAR_MASS/PAR_POSX/etc][ ParSortID[...] ]
//                3. Particles may be outside the target region
//                4. To ensure the consistency of random numbers, one must call the random number generator for
//                   ALL particles, including those too far away to affect the target region
//                5. No need to worry about the periodic boundary condition here
//                   --> Particle positions have been remapped in FB_AdvanceDt()
//                6. CoarseFine[] records the coarse-fine boundaries along the 26 sibling directions, defined as
//                            24  13  25
//                            15  05  17     z+1 plane
//                            22  12  23
//
//                            08  03  09
//                            00  XX  01     z   plane
//                            06  02  07
//                   y
//                   ^        20  11  21
//                   |        14  04  16     z-1 plane
//                   --->x    18  10  19
//                7. Invoked by FB_AdvanceDt()
//                8. Must NOT change particle positions
//                9. Since Fluid[] stores both the input and output data, the order of particles may affect the
//                   final output results
//                   --> For example, particle 2 may use the data updated by particle 1 as the input data
//                   --> Actually, even if we separate Fluid[] to input and output arrays, the final output results
//                       may still depend on the order of particles for non-local feedback since different particles
//                       may update the same cell
//                10. In general, it is recommended to have the maximum feedback radius no larger than half of the patch size
//                    (i.e., PATCH_SIZE/2=4 cells for PATCH_SIZE=8)
//                    --> Increase PATCH_SIZE if necessary
//                11. Linked to FB_User_Ptr in FB_Init_Plummer()
//                12. Must have FB_GHOST_SIZE>=1 for the mass accretion feedback
//
// Parameter   :  lv         : Target refinement level
//                TimeNew    : Target physical time to reach
//                TimeOld    : Physical time before update
//                             --> This function updates physical time from TimeOld to TimeNew
//                dt         : Time interval to advance solution
//                NPar       : Number of particles
//                ParSortID  : Sorted particle IDs
//                ParAtt     : Particle attribute arrays
//                Fluid      : Array to store the input/output fluid data
//                             --> Array size is fixed to (FB_NXT)^3=(PS2+2*FB_GHOST_SIZE)^3
//                EdgeL      : Left edge of Fluid[]
//                             --> Right edge is given by EdgeL[]+FB_NXT*dh
//                dh         : Cell size of Fluid[]
//                CoarseFine : Coarse-fine boundaries along the 26 sibling directions
//                TID        : Thread ID
//                RNG        : Random number generator
//                             --> Random number can be obtained by "RNG->GetValue( TID, Min, Max )",
//                                 where Min/Max specify the range of random numbers
//
// Return      :  Fluid, ParAtt
//-------------------------------------------------------------------------------------------------------
int FB_SinkAccretion( const int lv, const double TimeNew, const double TimeOld, const double dt,
                const int NPar, const int *ParSortID, real *ParAtt[PAR_NATT_TOTAL],
                real (*Fluid)[FB_NXT][FB_NXT][FB_NXT], const double EdgeL[], const double dh, bool CoarseFine[],
                const int TID, RandomNumber_t *RNG )
{

// check
#  ifdef GAMER_DEBUG
   if ( Fluid == NULL )    Aux_Error( ERROR_INFO, "Fluid == NULL !!\n" );
   if ( NPar > 0 )
   {
      if ( ParSortID == NULL )   Aux_Error( ERROR_INFO, "ParSortID == NULL for NPar = %d !!\n", NPar );
      if ( ParAtt == NULL )      Aux_Error( ERROR_INFO, "ParAtt == NULL for NPar = %d !!\n", NPar );
   }
#  endif // #ifdef GAMER_DEBUG

#  ifdef MY_DEBUG
   const char  FileName[] = "Record__Par_Acc";
   FILE *File = fopen( FileName, "a" );
#  endif

   real GasDens, DeltaM, Eg, Eg2, Ekin, Cell2Sinki, Cell2Sinkj, Cell2Sinkk, Cell2Sink2, GasRelVel[3]; 
   real ControlPosi[3], ControlPosj[3], ControlPosk[3];
   real Corner_Array[3]; // the corner of the ghost zone
   real GasMFracLeft;

   const int    AccCellNum     = 4;
   const int    MaxRemovalGas  = 100000;

   const double AccRadius      = AccCellNum*dh;
   const double _dh     = 1.0 / dh;
   const double dv      = CUBE( dh );
   const double _dv     = 1.0 / dv;
   const double epsilon = 0.001*dh;

   int      NRemoval           = 0;
   real   (*GasMFracLeftArr)   = new real [MaxRemovalGas];
   long   (*GasRemovalIdx)[3]  = new long [MaxRemovalGas][3];
   real   (*ParGain)[4]        = new real [NPar][4];

// prepare the corner array
   for (int d=0; d<3; d++)    Corner_Array[d] = EdgeL[d] + 0.5*dh ;

   for (int t=0; t<NPar; t++)
   {
      const int    p      = ParSortID[t];
      const double xyz[3] = { ParAtt[PAR_POSX][p], ParAtt[PAR_POSY][p], ParAtt[PAR_POSZ][p] }; // particle position
      real DeltaMSum      = 0;
      real DeltaMomSum[3] = { (real)0.0, (real)0.0, (real)0.0};

      int idx[3]; // cell idx in FB_NXT^3
      for (int d=0; d<3; d++)    idx[d] = (int)floor( ( xyz[d] - EdgeL[d] )*_dh );

      if ( idx[0] < FB_GHOST_SIZE-AccCellNum  ||  idx[0] >= FB_GHOST_SIZE+PS2+AccCellNum  ||
           idx[1] < FB_GHOST_SIZE-AccCellNum  ||  idx[1] >= FB_GHOST_SIZE+PS2+AccCellNum  ||
           idx[2] < FB_GHOST_SIZE-AccCellNum  ||  idx[2] >= FB_GHOST_SIZE+PS2+AccCellNum   ) // we want completed control volume
         continue;

//    Check the control volume
      for (int vki=idx[2]-AccCellNum; vki<=idx[2]+AccCellNum; vki++)
      for (int vji=idx[1]-AccCellNum; vji<=idx[1]+AccCellNum; vji++)
      for (int vii=idx[0]-AccCellNum; vii<=idx[0]+AccCellNum; vii++) // loop the nearby cells, to find the cells inside the control volumne (v)
      {
//       Inside the accretion radius
//       ===========================================================================================================
         ControlPosi[0] = Corner_Array[0] + vii*dh;
         ControlPosi[1] = Corner_Array[1] + vji*dh;
         ControlPosi[2] = Corner_Array[2] + vki*dh;
         Cell2Sinki = SQRT(SQR(ControlPosi[0] - xyz[0])+SQR(ControlPosi[1] - xyz[1])+SQR(ControlPosi[2] - xyz[2])); // distance to the sink
         if ( Cell2Sinki > AccRadius )                 continue; // check whether it is inside the control volume

//       Density threshold
//       ===========================================================================================================
         GasDens = Fluid[DENS][vki][vji][vii];
         if ( GasDens <= AccGasDensThres )                continue;

         DeltaM = (GasDens - AccGasDensThres)*dv; // the mass to be accreted

//       Central cell check
//       ===========================================================================================================
         bool NotCentralCell = true;
         if (( idx[0] == vii ) && ( idx[1] == vji ) && ( idx[2] == vki ))   
         NotCentralCell = false; // if pass, the following checks are skipped

         if ( NotCentralCell )
         {
//       Negative radial velocity
//       ===========================================================================================================
            GasRelVel[0] = Fluid[MOMX][vki][vji][vii]/GasDens - ParAtt[PAR_VELX][p];
            GasRelVel[1] = Fluid[MOMY][vki][vji][vii]/GasDens - ParAtt[PAR_VELY][p];
            GasRelVel[2] = Fluid[MOMZ][vki][vji][vii]/GasDens - ParAtt[PAR_VELZ][p];

            if ( GasRelVel[0] >= 0 ||  GasRelVel[1] >= 0 || GasRelVel[2] >= 0 )
            continue;

//       Bound state check
//       ===========================================================================================================
            real SelfPhi = (real)0.0; // self-potential
            for (int vkj=idx[2]-AccCellNum; vkj<=idx[2]+AccCellNum; vkj++)
            for (int vjj=idx[1]-AccCellNum; vjj<=idx[1]+AccCellNum; vjj++)
            for (int vij=idx[0]-AccCellNum; vij<=idx[0]+AccCellNum; vij++) // loop the nearby cells, to find the cells inside the control volumne (v)
            {  
               ControlPosj[0] = Corner_Array[0] + vij*dh;
               ControlPosj[1] = Corner_Array[1] + vjj*dh;
               ControlPosj[2] = Corner_Array[2] + vkj*dh;

               real rij = SQRT(SQR(ControlPosj[0] - ControlPosi[0])+SQR(ControlPosj[1] - ControlPosi[1])+SQR(ControlPosj[2] - ControlPosi[2]));
               if ( rij == 0.0 )                        continue;

               Cell2Sinkj = SQRT(SQR(ControlPosj[0] - xyz[0])+SQR(ControlPosj[1] - xyz[1])+SQR(ControlPosj[2] - xyz[2])); // distance to the center cell
               if ( Cell2Sinkj > AccRadius )            continue; // check whether it is inside the control volume

               SelfPhi += -NEWTON_G*Fluid[DENS][vkj][vjj][vij]*dv/rij; // potential
            } // vij, vjj, vkj

            SelfPhi += -NEWTON_G*ParAtt[PAR_MASS][p]/(Cell2Sinki+epsilon); // potential from the sink

            Eg   = DeltaM*SelfPhi; // no double counting here, since i is fixed
            Ekin = 0.5*DeltaM*( SQR(GasRelVel[0]) + SQR(GasRelVel[1]) + SQR(GasRelVel[2]));

            if (( Eg + Ekin ) >= 0)                     continue;

//       Overlapped accretion radius check
//       ===========================================================================================================
            bool NotMinEg = false;
            for (int tt=0; tt<NPar; tt++) // find the nearby sink
            {
               const int    pp      = ParSortID[tt];
               if ( pp == p )              continue;

               const double xxyyzz[3] = { ParAtt[PAR_POSX][pp], ParAtt[PAR_POSY][pp], ParAtt[PAR_POSZ][pp] }; // particle position
               Cell2Sink2 = SQRT(SQR(ControlPosi[0] - xxyyzz[0])+SQR(ControlPosi[1] - xxyyzz[1])+SQR(ControlPosi[2] - xxyyzz[2])); // distance to the sink
               if ( Cell2Sink2 > AccRadius )       continue;
               
               int idxx[3]; // cell idx in FB_NXT^3
               for (int d=0; d<3; d++)    idxx[d] = (int)floor( ( xxyyzz[d] - EdgeL[d] )*_dh );

               if ( idxx[0] < FB_GHOST_SIZE-AccCellNum  ||  idxx[0] >= FB_GHOST_SIZE+PS2+AccCellNum  ||
                     idxx[1] < FB_GHOST_SIZE-AccCellNum  ||  idxx[1] >= FB_GHOST_SIZE+PS2+AccCellNum  ||
                     idxx[2] < FB_GHOST_SIZE-AccCellNum  ||  idxx[2] >= FB_GHOST_SIZE+PS2+AccCellNum   ) // we want completed control volume
                     continue;

               real SelfPhi2 = (real)0.0; // self-potential
               for (int vkk=idxx[2]-AccCellNum; vkk<=idxx[2]+AccCellNum; vkk++)
               for (int vjk=idxx[1]-AccCellNum; vjk<=idxx[1]+AccCellNum; vjk++)
               for (int vik=idxx[0]-AccCellNum; vik<=idxx[0]+AccCellNum; vik++) // loop the nearby cells, to find the cells inside the control volumne (v)
               {
                  ControlPosk[0] = Corner_Array[0] + vik*dh;
                  ControlPosk[1] = Corner_Array[1] + vjk*dh;
                  ControlPosk[2] = Corner_Array[2] + vkk*dh;

                  real rik = SQRT(SQR(ControlPosk[0] - ControlPosi[0])+SQR(ControlPosk[1] - ControlPosi[1])+SQR(ControlPosk[2] - ControlPosi[2]));
                  if ( rik == 0.0 )                        continue;

                  Cell2Sinkk = SQRT(SQR(ControlPosk[0] - xxyyzz[0])+SQR(ControlPosk[1] - xxyyzz[1])+SQR(ControlPosk[2] - xxyyzz[2])); // distance to the center cell
                  if ( Cell2Sinkk > AccRadius )            continue; // check whether it is inside the control volume

                  SelfPhi2 += -NEWTON_G*Fluid[DENS][vkk][vjk][vik]*dv/rik; // potential
               } // vik, vjk, vkk

               SelfPhi2 += -NEWTON_G*ParAtt[PAR_MASS][pp]/(Cell2Sink2+epsilon); // potential from the sink

               Eg2   = DeltaM*SelfPhi2;
               if ( Eg2 < Eg )
               {
                  NotMinEg = true;
                  break;
               }
            } // for (int tt=0; tt<NPar; tt++)

            if ( NotMinEg )                              continue; 

         } // if ( CentralCell == false )

         DeltaMSum += DeltaM;
         GasMFracLeft = AccGasDensThres/GasDens;

         DeltaMomSum[0] += (1.0 - GasMFracLeft)*Fluid[MOMX][vki][vji][vii]*dv; // the momentum density of DeltaM
         DeltaMomSum[1] += (1.0 - GasMFracLeft)*Fluid[MOMY][vki][vji][vii]*dv;
         DeltaMomSum[2] += (1.0 - GasMFracLeft)*Fluid[MOMZ][vki][vji][vii]*dv;

         GasMFracLeftArr[NRemoval]  = GasMFracLeft;
         GasRemovalIdx[NRemoval][0] = vii;
         GasRemovalIdx[NRemoval][1] = vji;
         GasRemovalIdx[NRemoval][2] = vki;

#  ifdef MY_DEBUG
         fprintf( File,"%d %d %d %d %5.7e %5.7e", CentralCell, vii, vji, vki, GasMFracLeft*Fluid[DENS][vki][vji][vii], DeltaM);
         fprintf( File, "\n" );
#  endif

         NRemoval ++;
      } // vii, vji, vki

      ParGain[t][0] = DeltaMSum;
      ParGain[t][1] = (DeltaMomSum[0] + ParAtt[PAR_MASS][p]*ParAtt[PAR_VELX][p])/(DeltaMSum + ParAtt[PAR_MASS][p]);  // COM velocity of the sink after accretion
      ParGain[t][2] = (DeltaMomSum[1] + ParAtt[PAR_MASS][p]*ParAtt[PAR_VELY][p])/(DeltaMSum + ParAtt[PAR_MASS][p]);
      ParGain[t][3] = (DeltaMomSum[2] + ParAtt[PAR_MASS][p]*ParAtt[PAR_VELZ][p])/(DeltaMSum + ParAtt[PAR_MASS][p]);
   } // for (int t=0; t<NPar; t++)

// Remove the gas from Fluid
   for (int r=0; r<NRemoval; r++)
   {
      for (int v=0; v<NCOMP_TOTAL; v++)
      Fluid[v][GasRemovalIdx[r][2]][GasRemovalIdx[r][1]][GasRemovalIdx[r][0]] *= GasMFracLeftArr[r];
   } // for (int r=0; r<NRemoval; r++)

//  Update particle mass and velocity
//  ===========================================================================================================
   for (int t=0; t<NPar; t++)
   {
      const int    p      =  ParSortID[t];
      const double xyz[3] = { ParAtt[PAR_POSX][p], ParAtt[PAR_POSY][p], ParAtt[PAR_POSZ][p] }; // particle position

      int idx[3]; // cell idx in FB_NXT^3
      for (int d=0; d<3; d++)    idx[d] = (int)floor( ( xyz[d] - EdgeL[d] )*_dh );

      if ( idx[0] < FB_GHOST_SIZE  ||  idx[0] >= FB_GHOST_SIZE+PS2  ||
           idx[1] < FB_GHOST_SIZE  ||  idx[1] >= FB_GHOST_SIZE+PS2  ||
           idx[2] < FB_GHOST_SIZE  ||  idx[2] >= FB_GHOST_SIZE+PS2   )
         continue;

      ParAtt[PAR_MASS][p] += ParGain[t][0];
      ParAtt[PAR_VELX][p] =  ParGain[t][1];
      ParAtt[PAR_VELY][p] =  ParGain[t][2];
      ParAtt[PAR_VELZ][p] =  ParGain[t][3];
   }

#  ifdef MY_DEBUG
   // if ( NPar >= 1)
   // {
   //    fprintf( File,"TimeNew = %5.3e, NPar = %d", TimeNew, NPar);
   //    fprintf( File, "\n" );
   // }
   fclose( File );
#  endif

   delete [] GasMFracLeftArr;
   delete [] GasRemovalIdx;
   delete [] ParGain;

   return GAMER_SUCCESS;

} // FUNCTION : FB_SinkAccretion



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_End_Plummer
// Description :  Free the resources used by the user-specified feedback
//
// Note        :  1. Invoked by FB_End()
//                2. Linked to FB_End_User_Ptr in FB_Init_Plummer()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void FB_End_SinkAccretion()
{


} // FUNCTION : FB_End_SinkAccretion



//-------------------------------------------------------------------------------------------------------
// Function    :  FB_Init_Plummer
// Description :  Initialize the user-specified feedback
//
// Note        :  1. Invoked by FB_Init()
//                   --> Enable it by linking to the function pointer "FB_Init_User_Ptr"
//                2. Set FB_User_Ptr and FB_End_User_Ptr
//
// Parameter   :  None
//
// Return      :  FB_User_Ptr and FB_End_User_Ptr
//-------------------------------------------------------------------------------------------------------
void FB_Init_SinkAccretion()
{

   FB_User_Ptr     = FB_SinkAccretion;
   FB_End_User_Ptr = FB_End_SinkAccretion;

} // FUNCTION : FB_Init_SinkAccretion



#endif // #ifdef FEEDBACK