#include "CUFLU.h"



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "CUAPI.h"
#include "CUFLU_Shared_FluUtility.cu"
#include "CUDA_ConstMemory.h"

#endif // #ifdef __CUDACC__


// local function prototypes
#ifndef __CUDACC__

void Src_SetAuxArray_Cooling( double [], int [] );
void Src_SetConstMemory_Cooling( const double AuxArray_Flt[], const int AuxArray_Int[],
                                       double *&DevPtr_Flt, int *&DevPtr_Int );
void Src_SetFunc_Cooling( SrcFunc_t & );
void Src_WorkBeforeMajorFunc_Cooling( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                            double AuxArray_Flt[], int AuxArray_Int[] );
void Src_End_Cooling();

#endif



/********************************************************
1. Template of a user-defined source term
   --> Enabled by the runtime option "SRC_USER"

2. This file is shared by both CPU and GPU

   CUSRC_Src_Cooling.cu -> CPU_Src_Cooling.cpp

3. Four steps are required to implement a source term

   I.   Set auxiliary arrays
   II.  Implement the source-term function
   III. [Optional] Add the work to be done every time
        before calling the major source-term function
   IV.  Set initialization functions

4. The source-term function must be thread-safe and
   not use any global variable
********************************************************/



// =======================
// I. Set auxiliary arrays
// =======================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetAuxArray_Cooling
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by Src_Init_Cooling()
//                2. AuxArray_Flt/Int[] have the size of SRC_NAUX_USER defined in Macro.h (default = 10)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_SetAuxArray_Cooling( double AuxArray_Flt[], int AuxArray_Int[] )
{

   
   AuxArray_Flt[0] = ISM_UNIT_T2;
   AuxArray_Flt[1] = UNIT_T;
   AuxArray_Flt[2] = 1e-2;
   AuxArray_Flt[3] = 1e50;
   AuxArray_Flt[4] = GAMMA - 1.0;
   AuxArray_Flt[5] = 1 / AuxArray_Flt[4];
   AuxArray_Flt[6] = Const_kB;
   AuxArray_Flt[7] = POW( 10.0, 0.1 );
   AuxArray_Flt[8] = 0.2;

} // FUNCTION : Src_SetAuxArray_Cooling
#endif // #ifndef __CUDACC__



// ======================================
// II. Implement the source-term function
// ======================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Cooling
// Description :  Major source-term function
//
// Note        :  1. Invoked by CPU/GPU_SrcSolver_IterateAllCells()
//                2. See Src_SetAuxArray_Cooling() for the values stored in AuxArray_Flt/Int[]
//                3. Shared by both CPU and GPU
//
// Parameter   :  fluid             : Fluid array storing both the input and updated values
//                                    --> Including both active and passive variables
//                B                 : Cell-centered magnetic field
//                SrcTerms          : Structure storing all source-term variables
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                x/y/z             : Target physical coordinates
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                EoS               : EoS object
//                AuxArray_*        : Auxiliary arrays (see the Note above)
//
// Return      :  fluid[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void Src_Cooling( real fluid[], const real B[],
                               const SrcTerms_t *SrcTerms, const real dt, const real dh,
                               const double x, const double y, const double z,
                               const double TimeNew, const double TimeOld,
                               const real MinDens, const real MinPres, const real MinEint,
                               const EoS_t *EoS, const double AuxArray_Flt[], const int AuxArray_Int[] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif

   const real UNIT_T2          = AuxArray_Flt[0];
   const real Time2CGS         = AuxArray_Flt[1];
   const real Tmin             = AuxArray_Flt[2];
   const real Tmax             = AuxArray_Flt[3];
   const real Gamma_m1         = AuxArray_Flt[4];
   const real _Gamma_m1        = AuxArray_Flt[5];
   const real kB               = AuxArray_Flt[6];
   const real vardt            = AuxArray_Flt[7];
   const real varrel           = AuxArray_Flt[8];
   const bool CheckMinEint_Yes = true;
   const int Mode = 0;
   real Alpha_CT;
   real Etot, Emag, Eint_Code, Delta_Eint_Code;
   real Dens_Code, Dens_CGS;
   real Temp_CGS, Temp_Old_CGS;
   real Temp_Init_CGS, Temp_Update_CGS, Delta_Temp_CGS;
   real Temp_Code, Delta_Temp_Code;
   real CoolingRate0, CoolingRate1, dCoolingRatedTemp;
   real Delta_Time, Max_Delta_Time;
   real dt_CGS, Time = 0.0;
   real In_Flt[2] = {0};
   real Out[1] = {0};
   real Epsilon = 1e-5;
   int iterations = 0;

#  ifdef MHD
   Emag  = (real)0.5*( SQR(B[MAGX]) + SQR(B[MAGY]) + SQR(B[MAGZ]) );
#  else
   Emag  = (real)0.0;
#  endif
   
   Dens_Code = fluid[DENS];
   Eint_Code  = Hydro_Con2Eint( Dens_Code, fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], CheckMinEint_Yes, MinEint, Emag );
   Temp_Code  = Gamma_m1 * Eint_Code / Dens_Code;

   if ( Temp_Code <= 0.0 )
   {
      Temp_Init_CGS = Temp_Code * 1.4 * UNIT_T2;
      Temp_Update_CGS = 50.0;
      Delta_Temp_CGS = Temp_Update_CGS - Temp_Init_CGS;
      Delta_Temp_Code = Delta_Temp_CGS / UNIT_T2;
      Delta_Eint_Code = Delta_Temp_Code * Dens_Code * _Gamma_m1;
      fluid[ENGY] = Hydro_ConEint2Etot( Dens_Code, fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eint_Code + Delta_Eint_Code, Emag );
      return;
   }

   Temp_CGS   = Temp_Code * UNIT_T2;
   Temp_CGS  = MIN( MAX( Temp_CGS, Tmin), Tmax );
   Temp_Init_CGS= Temp_CGS;
   dt_CGS = dt * Time2CGS;
   
   if ( Dens_CGS <= 1e-10 )   Dens_CGS = 1e-10;
   Alpha_CT = Dens_CGS * kB * _Gamma_m1;

   while( Time < dt_CGS )
   {
      if ( Temp_CGS < 0.0 )
      {
         Temp_CGS = MIN( 4000 / MAX(Dens_CGS, 1e-10), 8000);
      }
      Temp_Old_CGS = Temp_CGS;
      In_Flt[0] = Dens_CGS;
      In_Flt[1] = Temp_CGS;
      EoS->General_FuncPtr(Mode, Out, In_Flt, NULL, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table);
      CoolingRate0 = Out[0];

      In_Flt[0] = Dens_CGS;
      In_Flt[1] = Temp_CGS * Epsilon;
      EoS->General_FuncPtr(Mode, Out, In_Flt, NULL, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table);
      CoolingRate1 = Out[0];

      dCoolingRatedTemp = (CoolingRate1 - CoolingRate0) / (Temp_CGS * Epsilon);
      
      if (iterations == 0)
      {
         if ( dCoolingRatedTemp != 0.0 )
         {
            Delta_Time = FABS( 0.1 * Alpha_CT / dCoolingRatedTemp );
         }
         else
         {
            Delta_Time = 0.1 * dt_CGS;
         }

         Max_Delta_Time = dt_CGS - Time;
         
         if ( Delta_Time > 0.7 * Max_Delta_Time )
         {
            Delta_Time = Max_Delta_Time * ( 1.0 + 1e-12 );
         }
      }

      Delta_Temp_CGS = CoolingRate0 / ( Alpha_CT / Delta_Time - dCoolingRatedTemp );
      Epsilon = FABS(Delta_Temp_CGS / Temp_CGS);
      if ( Epsilon > 0.2 )
      {
         Delta_Temp_CGS = 0.2 * Temp_Code * Delta_Temp_CGS / FABS(Delta_Temp_CGS);
      }

      iterations ++;
      Temp_CGS += Delta_Temp_CGS;
      Time += Delta_Time;

      Delta_Time = vardt * varrel * Delta_Time / MAX( vardt * Epsilon, varrel );
      Max_Delta_Time = dt_CGS - Time;   
      if ( Delta_Time > 0.7 * Max_Delta_Time )
      {
         Delta_Time = Max_Delta_Time * ( 1.0 + 1e-12 );
      }

      if ( Temp_CGS < 0.0 )
      {
         Temp_CGS = 100.0;
      }
   }
   Temp_Update_CGS = Temp_CGS;
   
   Delta_Temp_CGS = Temp_Update_CGS - Temp_Init_CGS;
   Delta_Temp_Code = Delta_Temp_CGS / UNIT_T2;
   Delta_Eint_Code = Delta_Temp_Code * Dens_Code * _Gamma_m1;
   fluid[ENGY] = Hydro_ConEint2Etot( Dens_Code, fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eint_Code + Delta_Eint_Code, Emag );

} // FUNCTION : Src_Cooling



// ==================================================
// III. [Optional] Add the work to be done every time
//      before calling the major source-term function
// ==================================================

//-------------------------------------------------------------------------------------------------------
// Function    :  Src_WorkBeforeMajorFunc_Cooling
// Description :  Specify work to be done every time before calling the major source-term function
//
// Note        :  1. Invoked by Src_WorkBeforeMajorFunc()
//                   --> By linking to "Src_WorkBeforeMajorFunc_User_Ptr" in Src_Init_Cooling()
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  lv               : Target refinement level
//                TimeNew          : Target physical time to reach
//                TimeOld          : Physical time before update
//                                   --> The major source-term function will update the system from TimeOld to TimeNew
//                dt               : Time interval to advance solution
//                                   --> Physical coordinates : TimeNew - TimeOld == dt
//                                       Comoving coordinates : TimeNew - TimeOld == delta(scale factor) != dt
//                AuxArray_Flt/Int : Auxiliary arrays
//                                   --> Can be used and/or modified here
//                                   --> Must call Src_SetConstMemory_Cooling() after modification
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void Src_WorkBeforeMajorFunc_Cooling( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                            double AuxArray_Flt[], int AuxArray_Int[] )
{

// uncomment the following lines if the auxiliary arrays have been modified
//#  ifdef GPU
//   Src_SetConstMemory_Cooling( AuxArray_Flt, AuxArray_Int,
//                                     SrcTerms.User_AuxArrayDevPtr_Flt, SrcTerms.User_AuxArrayDevPtr_Int );
//#  endif

} // FUNCTION : Src_WorkBeforeMajorFunc_Cooling
#endif



// ================================
// IV. Set initialization functions
// ================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE SrcFunc_t SrcFunc_Ptr = Src_Cooling;

//-----------------------------------------------------------------------------------------
// Function    :  Src_SetFunc_Cooling
// Description :  Return the function pointer of the CPU/GPU source-term function
//
// Note        :  1. Invoked by Src_Init_Cooling()
//                2. Call-by-reference
//                3. Use either CPU or GPU but not both of them
//
// Parameter   :  SrcFunc_CPU/GPUPtr : CPU/GPU function pointer to be set
//
// Return      :  SrcFunc_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void Src_SetFunc_Cooling( SrcFunc_t &SrcFunc_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &SrcFunc_GPUPtr, SrcFunc_Ptr, sizeof(SrcFunc_t) )  );
}

#elif ( !defined GPU )

void Src_SetFunc_Cooling( SrcFunc_t &SrcFunc_CPUPtr )
{
   SrcFunc_CPUPtr = SrcFunc_Ptr;
}

#endif // #ifdef __CUDACC__ ... elif ...



#ifdef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Src_SetConstMemory_Cooling
// Description :  Set the constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Src_Init_Cooling() and, if necessary, Src_WorkBeforeMajorFunc_Cooling()
//                3. SRC_NAUX_USER is defined in Macro.h
//
// Parameter   :  AuxArray_Flt/Int : Auxiliary arrays to be copied to the constant memory
//                DevPtr_Flt/Int   : Pointers to store the addresses of constant memory arrays
//
// Return      :  c_Src_User_AuxArray_Flt[], c_Src_User_AuxArray_Int[], DevPtr_Flt, DevPtr_Int
//---------------------------------------------------------------------------------------------------
void Src_SetConstMemory_Cooling( const double AuxArray_Flt[], const int AuxArray_Int[],
                                       double *&DevPtr_Flt, int *&DevPtr_Int )
{

// copy data to constant memory
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_Cooling_AuxArray_Flt, AuxArray_Flt, SRC_NAUX_COOLING*sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Src_Cooling_AuxArray_Int, AuxArray_Int, SRC_NAUX_COOLING*sizeof(int   ) )  );

// obtain the constant-memory pointers
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Flt, c_Src_Cooling_AuxArray_Flt) );
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&DevPtr_Int, c_Src_Cooling_AuxArray_Int) );

} // FUNCTION : Src_SetConstMemory_Cooling
#endif // #ifdef __CUDACC__



#ifndef __CUDACC__

//-----------------------------------------------------------------------------------------
// Function    :  Src_Init_Cooling
// Description :  Initialize a user-specified source term
//
// Note        :  1. Set auxiliary arrays by invoking Src_SetAuxArray_*()
//                   --> Copy to the GPU constant memory and store the associated addresses
//                2. Set the source-term function by invoking Src_SetFunc_*()
//                   --> Unlike other modules (e.g., EoS), here we use either CPU or GPU but not
//                       both of them
//                3. Set the function pointers "Src_WorkBeforeMajorFunc_User_Ptr" and "Src_End_User_Ptr"
//                4. Invoked by Src_Init()
//                   --> Enable it by linking to the function pointer "Src_Init_User_Ptr"
//                5. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_Init_Cooling()
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );
// set the auxiliary arrays
   Src_SetAuxArray_Cooling( Src_Cooling_AuxArray_Flt, Src_Cooling_AuxArray_Int );

// copy the auxiliary arrays to the GPU constant memory and store the associated addresses
#  ifdef GPU
   Src_SetConstMemory_Cooling( Src_Cooling_AuxArray_Flt, Src_Cooling_AuxArray_Int,
                               SrcTerms.Cooling_AuxArrayDevPtr_Flt, SrcTerms.Cooling_AuxArrayDevPtr_Int );
#  else
   SrcTerms.Cooling_AuxArrayDevPtr_Flt = Src_Cooling_AuxArray_Flt;
   SrcTerms.Cooling_AuxArrayDevPtr_Int = Src_Cooling_AuxArray_Int;
#  endif

// set the major source-term function
   Src_SetFunc_Cooling( SrcTerms.Cooling_FuncPtr );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
} // FUNCTION : Src_Init_Cooling



//-----------------------------------------------------------------------------------------
// Function    :  Src_End_Cooling
// Description :  Free the resources used by a user-specified source term
//
// Note        :  1. Invoked by Src_End()
//                   --> Enable it by linking to the function pointer "Src_End_User_Ptr"
//                2. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Src_End_Cooling()
{


} // FUNCTION : Src_End_Cooling

#endif // #ifndef __CUDACC__
