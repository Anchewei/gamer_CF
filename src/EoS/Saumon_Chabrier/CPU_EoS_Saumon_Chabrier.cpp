#include "EoS_Saumon_Chabrier.h"
#include "CUFLU.h"
#ifdef __CUDACC__
#include "CUAPI.h"
#include "CUFLU_Shared_FluUtility.cu"
#endif

#if ( MODEL == HYDRO && EOS == EOS_SAUMON_CHABRIER )

// global variables
int nRho;
int nEner;
int nTemp;
int nPres;
double MaxRho;
double MinRho;
double MaxEner;
double MinEner;
double MaxTemp;
double MinTemp;
double MaxPres;
double MinPres;
double dRho;
double dEner;
double dTemp;
double dPres;
double *Table_LogRho = NULL;
double *Table_LogEner = NULL;
double *Table_DensEint2Temp = NULL;
double *Table_DensEint2Pres = NULL;
double *Table_DensEint2Cs = NULL;
double *Table_DensPres2Eint = NULL;
double *Table_DensTemp2Eint = NULL;
bool flag = false;

// prototypes
void Read_Saumon_Chabrier_Table( char *ISM_Saumon_Chabrier_Table_Name );
void Fill_Hole( double *Table, const int Table_Index );
void Create_Table( const int Target_Table_Index );

/********************************************************
1. Template of a user-defined EoS (EOS_USER)

2. This file is shared by both CPU and GPU

   GPU_EoS_Saumon_Chabrier.cu -> CPU_EoS_Saumon_Chabrier.cpp

3. Three steps are required to implement an EoS

   I.   Set EoS auxiliary arrays
   II.  Implement EoS conversion functions
   III. Set EoS initialization functions

4. All EoS conversion functions must be thread-safe and
   not use any global variable

5. When an EoS conversion function fails, it is recommended
   to return NAN in order to trigger auto-correction such as
   "OPT__1ST_FLUX_CORR" and "AUTO_REDUCE_DT"
********************************************************/



// =============================================
// I. Set EoS auxiliary arrays
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_SetAuxArray_Saumon_Chabrier
// Description :  Set the auxiliary arrays AuxArray_Flt/Int[]
//
// Note        :  1. Invoked by EoS_Init_Saumon_Chabrier()
//                2. AuxArray_Flt/Int[] have the size of EOS_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//                4. Physical constants such as Const_amu/Const_kB should be set to unity when disabling OPT__UNIT
//
// Parameter   :  AuxArray_Flt/Int : Floating-point/Integer arrays to be filled up
//
// Return      :  AuxArray_Flt/Int[]
//-------------------------------------------------------------------------------------------------------
#ifndef __CUDACC__
void EoS_SetAuxArray_Saumon_Chabrier( double AuxArray_Flt[], int AuxArray_Int[] )
{

   AuxArray_Flt[SAUMON_CHABRIER_AUX_MAXRHO] = MaxRho;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_MINRHO] = MinRho;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_MAXENER] = MaxEner;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_MINENER] = MinEner;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_MAXTEMP] = MaxTemp;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_MINTEMP] = MinTemp;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_MAXPRES] = MaxPres;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_MINPRES] = MinPres;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_DRHO] = dRho;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_DENER] = dEner;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_DTEMP] = dTemp;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_DPRES] = dPres;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_D] = UNIT_D;
   AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_V] = UNIT_V;
   
   AuxArray_Int[SAUMON_CHABRIER_AUX_NRHO] = nRho;
   AuxArray_Int[SAUMON_CHABRIER_AUX_NENER] = nEner;
   AuxArray_Int[SAUMON_CHABRIER_AUX_NTEMP] = nTemp;
   AuxArray_Int[SAUMON_CHABRIER_AUX_NPRES] = nPres;

} // FUNCTION : EoS_SetAuxArray_Saumon_Chabrier
#endif // #ifndef __CUDACC__



// =============================================
// II. Implement EoS conversion functions
//     (1) EoS_DensEint2Pres_*
//     (2) EoS_DensPres2Eint_*
//     (3) EoS_DensPres2CSqr_*
//     (4) EoS_DensEint2Temp_* [OPTIONAL]
//     (5) EoS_DensTemp2Pres_* [OPTIONAL]
//     (6) EoS_DensEint2Entr_* [OPTIONAL]
//     (7) EoS_General_*       [OPTIONAL]
// =============================================

//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Pres_Saumon_Chabrier
// Description :  Convert gas mass density and internal energy density to gas pressure
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_Saumon_Chabrier() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Pres_Saumon_Chabrier( const real Dens, const real Eint, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Dens, "input density",         ERROR_INFO, UNPHY_VERBOSE );
// note that some EoS may support Eint<0
   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Eint, "input internal energy", ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG
   
   if ( Eint <= 0.0 || Dens <= 0.0 ) return 0.0;

   const real UNITD = (real) AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_D];
   const real UNITV = (real) AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_V];
   const real d_Rho = (real) AuxArray_Flt[SAUMON_CHABRIER_AUX_DRHO];
   const real d_Ener = (real) AuxArray_Flt[SAUMON_CHABRIER_AUX_DENER];
   const real Min_Rho = (real) AuxArray_Flt[SAUMON_CHABRIER_AUX_MINRHO];
   const real Min_Ener = (real) AuxArray_Flt[SAUMON_CHABRIER_AUX_MINENER];
   const int NRHO = AuxArray_Int[SAUMON_CHABRIER_AUX_NRHO];
   const int NENER = AuxArray_Int[SAUMON_CHABRIER_AUX_NENER];
   
   const real Dens_CGS = Dens * UNITD;
   const real Eint_CGS = Eint * UNITD * SQR( UNITV );
   const real lr = (LOG10(Dens_CGS) - Min_Rho) / d_Rho;
   const real le = (LOG10(Eint_CGS) - Min_Ener - LOG10(Dens_CGS)) / d_Ener;
   if ( lr < 0 || lr > NRHO - 1 ) return 0.0; // Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e, MinRho = %20.14e, dRho = %20.14e\n", Dens_CGS, Min_Rho, d_Rho );
   if ( le < 0 || le > NENER - 1 ) return 0.0; //Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e, Eint = %20.14e, MinEner = %20.14e, dEner = %20.14e\n", Dens_CGS, Eint_CGS, Min_Ener, d_Ener );
   
   const int ir = FLOOR(lr);
   const int ie = FLOOR(le);
   const real x = lr - ir;
   const real y = le - ie;

   real Pres_CGS = 0.0;
   real Pres_Code = 0.0;
   Pres_CGS += ( 1 - x ) * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ ie * NRHO + ir ];
   Pres_CGS +=       x   * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ ie * NRHO + ir + 1 ];
   Pres_CGS += ( 1 - x ) *       y   * Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ (ie + 1)* NRHO + ir ];
   Pres_CGS +=       x   *       y   * Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ (ie + 1)* NRHO + ir + 1];
   Pres_Code = Pres_CGS / (UNITD * SQR(UNITV) );

   if ( Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ ie * NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ ie * NRHO + ir + 1 ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ (ie + 1)* NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ (ie + 1)* NRHO + ir + 1] == 0.0 )
   {
      Pres_Code = 0.0;
   }

// check
#  ifdef GAMER_DEBUG
   if ( Hydro_CheckUnphysical( UNPHY_MODE_SING, &Pres, "output pressure", ERROR_INFO, UNPHY_VERBOSE ) )
   {
      printf( "Dens=%13.7e, Eint=%13.7e\n", Dens, Eint );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Pres_Code;

} // FUNCTION : EoS_DensEint2Pres_Saumon_Chabrier



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2Eint_Saumon_Chabrier
// Description :  Convert gas mass density and pressure to gas internal energy density
//
// Note        :  1. See EoS_DensEint2Pres_Saumon_Chabrier()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas internal energy density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2Eint_Saumon_Chabrier( const real Dens, const real Pres, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Dens, "input density",  ERROR_INFO, UNPHY_VERBOSE );
   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Pres, "input pressure", ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG

   if ( Dens <= 0.0 || Pres <= 0.0 ) return 0.0;

   const real UNITD = AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_D];
   const real UNITV = AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_V];
   const real d_Rho = AuxArray_Flt[SAUMON_CHABRIER_AUX_DRHO];
   const real d_Pres = AuxArray_Flt[SAUMON_CHABRIER_AUX_DPRES];
   const real Min_Rho = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINRHO];
   const real Min_Pres = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINPRES];
   const int NRHO = AuxArray_Int[SAUMON_CHABRIER_AUX_NRHO];
   const int NPRES = AuxArray_Int[SAUMON_CHABRIER_AUX_NPRES];
   const real Dens_CGS = Dens * UNITD;
   const real Pres_CGS = Pres * UNITD * SQR(UNITV);
   const real lr = ( LOG10(Dens_CGS) - Min_Rho ) / d_Rho;
   const real lp = ( LOG10(Pres_CGS) - LOG10(Min_Pres) - LOG10(Dens_CGS) ) / d_Pres;
   if ( lr < 0 || lr > NRHO - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e, MinRho = %20.14e, dRho = %20.14e\n", Dens_CGS, Min_Rho, d_Rho );
   if ( lp < 0 || lp > NPRES - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e, Pres = %20.14e, MinPres = %20.14e, dPres = %20.14e\n", Dens_CGS, Pres_CGS, Min_Pres, d_Pres );
   
   const int ir = FLOOR(lr);
   const int ip = FLOOR(lp);
   const real x = lr - ir;
   const real y = lp - ip;

   real Eint_CGS = 0.0;
   real Eint_Code = 0.0;
   Eint_CGS += ( 1 - x ) * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ ip * NRHO + ir ];
   Eint_CGS +=       x   * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ ip * NRHO + ir + 1 ];
   Eint_CGS += ( 1 - x ) *       y   * Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ (ip + 1)* NRHO + ir ];
   Eint_CGS +=       x   *       y   * Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ (ip + 1)* NRHO + ir + 1];
   Eint_Code = Eint_CGS / (UNITD * SQR(UNITV) );

   if ( Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ ip * NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ ip * NRHO + ir + 1 ] *
        Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ (ip + 1)* NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ (ip + 1)* NRHO + ir + 1] == 0.0 )
   {
      Eint_Code = 0.0;
   }


// check
#  ifdef GAMER_DEBUG
// note that some EoS may support Eint<0
   if ( Hydro_CheckUnphysical( UNPHY_MODE_SING, &Eint, "output internal energy", ERROR_INFO, UNPHY_VERBOSE ) )
   {
      printf( "Dens=%13.7e, Pres=%13.7e\n", Dens, Pres );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Eint_Code;

} // FUNCTION : EoS_DensPres2Eint_Saumon_Chabrier



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensPres2CSqr_Saumon_Chabrier
// Description :  Convert gas mass density and pressure to sound speed squared
//
// Note        :  1. See EoS_DensEint2Pres_Saumon_Chabrier()
//
// Parameter   :  Dens       : Gas mass density
//                Pres       : Gas pressure
//                Passive    : Passive scalars
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Sound speed squared
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensPres2CSqr_Saumon_Chabrier( const real Dens, const real Pres, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  if ( NCOMP_PASSIVE > 0 )
   if ( Passive == NULL )  printf( "ERROR : Passive == NULL in %s !!\n", __FUNCTION__ );
#  endif
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Dens, "input density" , ERROR_INFO, UNPHY_VERBOSE );
   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Pres, "input pressure", ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG
   if ( Dens <= 0.0 || Pres <= 0.0 ) return 0.0;

   const real UNITD = AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_D];
   const real UNITV = AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_V];
   const real d_Rho = AuxArray_Flt[SAUMON_CHABRIER_AUX_DRHO];
   const real d_Pres = AuxArray_Flt[SAUMON_CHABRIER_AUX_DPRES];
   const real d_Ener = AuxArray_Flt[SAUMON_CHABRIER_AUX_DENER];
   const real Min_Rho = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINRHO];
   const real Min_Pres = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINPRES];
   const real Min_Ener = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINENER];
   const int NRHO = AuxArray_Int[SAUMON_CHABRIER_AUX_NRHO];
   const int NENER = AuxArray_Int[SAUMON_CHABRIER_AUX_NENER];
   const int NPRES = AuxArray_Int[SAUMON_CHABRIER_AUX_NPRES];
   const real Dens_CGS = Dens * UNITD;
   const real Pres_CGS = Pres * UNITD * SQR(UNITV);
   const real lr = ( LOG10(Dens_CGS) - Min_Rho) / d_Rho;
   const real lp = ( LOG10(Pres_CGS) - LOG10(Min_Pres) - LOG10(Dens_CGS) ) / d_Pres;
   if ( lr < 0 || lr > NRHO - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e, MinRho = %20.14e, dRho = %20.14e\n", Dens_CGS, Min_Rho, d_Rho );
   if ( lp < 0 || lp > NPRES - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e, Pres = %20.14e, MinPres = %20.14e, dPres = %20.14e\n", Dens_CGS, Pres_CGS, Min_Pres, d_Pres );
   
   const int ir = FLOOR(lr);
   const int ip = FLOOR(lp);
   real x = lr - ir;
   real y = lp - ip;

   real Eint_CGS = 0.0;
   real Eint_Code = 0.0;
   real Cs_CGS = 0.0;
   real Cs2_Code = 0.0;

   Eint_CGS += ( 1 - x ) * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ ip * NRHO + ir ];
   Eint_CGS +=       x   * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ ip * NRHO + ir + 1 ];
   Eint_CGS += ( 1 - x ) *       y   * Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ (ip + 1)* NRHO + ir ];
   Eint_CGS +=       x   *       y   * Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ (ip + 1)* NRHO + ir + 1];
   Eint_Code =  Eint_CGS / ( UNITD * SQR(UNITV) );

   if ( Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ ip * NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ ip * NRHO + ir + 1 ] *
        Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ (ip + 1)* NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT][ (ip + 1)* NRHO + ir + 1] == 0.0 )
   {
      Eint_Code = 0.0;
   }
   
   Eint_CGS = Eint_Code * UNITD * SQR(UNITV);
   if ( Eint_CGS <= 0.0 ) return 0.0;
   
   const real le = ( LOG10(Eint_CGS) - Min_Ener - LOG10(Dens_CGS) ) / d_Ener;
   if ( le < 0 || le > NENER - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e, Eint = %20.14e, MinEner = %20.14e, dEner = %20.14e\n", Dens_CGS, Eint_CGS, Min_Ener, d_Ener );   
   
   const int ie = FLOOR(le);
   x = lr - ir;
   y = le - ie;
   
   Cs_CGS += ( 1 - x ) * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSEINT2CS][ ie * NRHO + ir ];
   Cs_CGS +=       x   * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSEINT2CS][ ie * NRHO + ir + 1 ];
   Cs_CGS += ( 1 - x ) *       y   * Table[SAUMON_CHABRIER_TAB_DENSEINT2CS][ (ie + 1)* NRHO + ir ];
   Cs_CGS +=       x   *       y   * Table[SAUMON_CHABRIER_TAB_DENSEINT2CS][ (ie + 1)* NRHO + ir + 1];
   Cs2_Code = SQR( Cs_CGS / UNIT_V );

   if ( Table[SAUMON_CHABRIER_TAB_DENSEINT2CS][ ie * NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2CS][ ie * NRHO + ir + 1 ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2CS][ (ie + 1)* NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2CS][ (ie + 1)* NRHO + ir + 1] == 0.0 )
   {
      Cs2_Code = 0.0;
   }

// check
#  ifdef GAMER_DEBUG
   if ( Hydro_CheckUnphysical( UNPHY_MODE_SING, &Cs2, "output sound speed squared", ERROR_INFO, UNPHY_VERBOSE ) )
   {
      printf( "Dens=%13.7e, Pres=%13.7e\n", Dens, Pres );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Cs2_Code;

} // FUNCTION : EoS_DensPres2CSqr_Saumon_Chabrier



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Temp_Saumon_Chabrier
// Description :  Convert gas mass density and internal energy density to gas temperature
//
// Note        :  1. Internal energy density here is per unit volume instead of per unit mass
//                2. See EoS_SetAuxArray_Saumon_Chabrier() for the values stored in AuxArray_Flt/Int[]
//                3. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas temperature in kelvin
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Temp_Saumon_Chabrier( const real Dens, const real Eint, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Dens, "input density"        , ERROR_INFO, UNPHY_VERBOSE );
   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Eint, "input internal energy", ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG

   if ( Dens <= 0.0 || Eint <= 0.0 ) return 0.0;
   const real UNITD = AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_D];
   const real UNITV = AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_V];
   const real d_Rho = AuxArray_Flt[SAUMON_CHABRIER_AUX_DRHO];
   const real d_Ener = AuxArray_Flt[SAUMON_CHABRIER_AUX_DENER];
   const real Min_Rho = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINRHO];
   const real Min_Ener = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINENER];
   const int NRHO = AuxArray_Int[SAUMON_CHABRIER_AUX_NRHO];
   const int NENER = AuxArray_Int[SAUMON_CHABRIER_AUX_NENER];
   const real Dens_CGS = Dens * UNITD;
   const real Eint_CGS = Eint * UNITD * SQR( UNITV );
   const real lr = (LOG10(Dens_CGS) - Min_Rho) / d_Rho;
   const real le = (LOG10(Eint_CGS) - Min_Ener - LOG10(Dens_CGS) ) / d_Ener;
   if ( lr < 0 || lr > NRHO - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e\n", Dens_CGS );
   if ( le < 0 || le > NENER - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Eint = %20.14e\n", Eint_CGS );
   
   const int ir = FLOOR(lr);
   const int ie = FLOOR(le);
   const real x = lr - ir;
   const real y = le - ie;

   real Temp = 0.0;
   Temp += ( 1.0 - x ) * ( 1.0 - y ) * Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP][ ie * NRHO + ir ];
   Temp +=       x     * ( 1.0 - y ) * Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP][ ie * NRHO + ir + 1 ];
   Temp += ( 1.0 - x ) *         y   * Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP][ (ie + 1)* NRHO + ir ];
   Temp +=       x     *         y   * Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP][ (ie + 1)* NRHO + ir + 1];

   if ( Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP][ ie * NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP][ ie * NRHO + ir + 1 ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP][ (ie + 1)* NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP][ (ie + 1)* NRHO + ir + 1] == 0.0 )
   {
      Temp = 0.0;
   }


// check
#  ifdef GAMER_DEBUG
   if ( Hydro_CheckUnphysical( UNPHY_MODE_SING, &Temp, "output temperature", ERROR_INFO, UNPHY_VERBOSE ) )
   {
      printf( "Dens=%13.7e, Eint=%13.7e\n", Dens, Eint );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Temp;

} // FUNCTION : EoS_DensEint2Temp_Saumon_Chabrier



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensTemp2Pres_Saumon_Chabrier
// Description :  Convert gas mass density and temperature to gas pressure
//
// Note        :  1. See EoS_SetAuxArray_Saumon_Chabrier() for the values stored in AuxArray_Flt/Int[]
//                2. Temperature is in kelvin
//
// Parameter   :  Dens       : Gas mass density
//                Temp       : Gas temperature in kelvin
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensTemp2Pres_Saumon_Chabrier( const real Dens, const real Temp, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Dens, "input density"    , ERROR_INFO, UNPHY_VERBOSE );
   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Temp, "input temperature", ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG

   const real UNITD = AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_D];
   const real UNITV = AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_V];
   const real d_Rho = AuxArray_Flt[SAUMON_CHABRIER_AUX_DRHO];
   const real d_Temp = AuxArray_Flt[SAUMON_CHABRIER_AUX_DTEMP];
   const real d_Ener = AuxArray_Flt[SAUMON_CHABRIER_AUX_DENER];
   const real Min_Rho = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINRHO];
   const real Min_Temp = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINTEMP];
   const real Min_Ener = AuxArray_Flt[SAUMON_CHABRIER_AUX_MINENER];
   const int NRHO = AuxArray_Int[SAUMON_CHABRIER_AUX_NRHO];
   const int NTEMP = AuxArray_Int[SAUMON_CHABRIER_AUX_NTEMP];
   const int NENER = AuxArray_Int[SAUMON_CHABRIER_AUX_NENER];
   const real Dens_CGS = Dens * UNITD;
   const real lr = ( LOG10(Dens_CGS) - Min_Rho ) / d_Rho;
   const real lt = ( LOG10(Temp) - LOG10(Min_Temp) ) / d_Temp;
   if ( lr < 0 || lr > NRHO - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Rho = %20.14e\n", Dens_CGS );
   if ( lt < 0 || lt > NTEMP - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Temp = %20.14e\n", Temp );
   
   const int ir = FLOOR(lr);
   const int it = FLOOR(lt);
   real x = lr - ir;
   real y = lt - it;
   
   real Eint_CGS = 0.0;
   real Eint_Code = 0.0;
   Eint_CGS += ( 1 - x ) * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT][ it * NRHO + ir ];
   Eint_CGS +=       x   * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT][ it * NRHO + ir + 1 ];
   Eint_CGS += ( 1 - x ) *       y   * Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT][ (it + 1)* NRHO + ir ];
   Eint_CGS +=       x   *       y   * Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT][ (it + 1)* NRHO + ir + 1];
   Eint_Code = Eint_CGS / ( UNITD * SQR(UNITV) );

   if ( Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT][ it * NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT][ it * NRHO + ir + 1 ] *
        Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT][ (it + 1)* NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT][ (it + 1)* NRHO + ir + 1] == 0.0 )
   {
      Eint_Code = 0.0;
   }

   Eint_CGS = Eint_Code * UNITD * SQR(UNITV);
   if ( Eint_CGS == 0.0 ) return 0.0;
   
   const real le = (LOG10(Eint_CGS) - Min_Ener - LOG10(Dens_CGS) ) / d_Ener;
   if ( le < 0 || le > NENER - 1 ) Aux_Error(ERROR_INFO, "Interpolation Out of Bound, Eint = %20.14e\n", Eint_CGS );
   
   const int ie = FLOOR(le);
   x = lr - ir;
   y = le - ie;

   real Pres_CGS = 0.0;
   real Pres_Code = 0.0;
   Pres_CGS += ( 1 - x ) * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ ie * NRHO + ir ];
   Pres_CGS +=       x   * ( 1 - y ) * Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ ie * NRHO + ir + 1 ];
   Pres_CGS += ( 1 - x ) *       y   * Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ (ie + 1)* NRHO + ir ];
   Pres_CGS +=       x   *       y   * Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ (ie + 1)* NRHO + ir + 1];
   Pres_Code = Pres_CGS / ( UNITD * SQR(UNITV) );

   if ( Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ ie * NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ ie * NRHO + ir + 1 ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ (ie + 1)* NRHO + ir ] *
        Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES][ (ie + 1)* NRHO + ir + 1] == 0.0 )
   {
      Pres_Code = 0.0;
   }


// check
#  ifdef GAMER_DEBUG
   if ( Hydro_CheckUnphysical( UNPHY_MODE_SING, &Pres, "output pressure", ERROR_INFO, UNPHY_VERBOSE ) )
   {
      printf( "Dens=%13.7e, Temp=%13.7e\n", Dens, Temp );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Pres_Code;

} // FUNCTION : EoS_DensTemp2Pres_Saumon_Chabrier



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_DensEint2Entr_Saumon_Chabrier
// Description :  Convert gas mass density and internal energy density to gas entropy
//
// Note        :  1. See EoS_SetAuxArray_Saumon_Chabrier() for the values stored in AuxArray_Flt/Int[]
//
// Parameter   :  Dens       : Gas mass density
//                Eint       : Gas internal energy density
//                Passive    : Passive scalars (must not used here)
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Gas entropy
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static real EoS_DensEint2Entr_Saumon_Chabrier( const real Dens, const real Eint, const real Passive[],
                                             const double AuxArray_Flt[], const int AuxArray_Int[],
                                             const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );

   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Dens, "input density"        , ERROR_INFO, UNPHY_VERBOSE );
   Hydro_CheckUnphysical( UNPHY_MODE_SING, &Eint, "input internal energy", ERROR_INFO, UNPHY_VERBOSE );
#  endif // GAMER_DEBUG


   real Entr = -1.0;

   /*
   Entr = ...;
   */


// check
#  ifdef GAMER_DEBUG
   if ( Hydro_CheckUnphysical( UNPHY_MODE_SING, &Entr, "output entropy", ERROR_INFO, UNPHY_VERBOSE ) )
   {
      printf( "Dens=%13.7e, Eint=%13.7e\n", Dens, Eint );
#     if ( NCOMP_PASSIVE > 0 )
      printf( "Passive scalars:" );
      for (int v=0; v<NCOMP_PASSIVE; v++)    printf( " %d=%13.7e", v, Passive[v] );
      printf( "\n" );
#     endif
   }
#  endif // GAMER_DEBUG


   return Entr;

} // FUNCTION : EoS_DensEint2Entr_Saumon_Chabrier



//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_General_Saumon_Chabrier
// Description :  General EoS converter: In_*[] -> Out[]
//
// Note        :  1. See EoS_DensEint2Pres_Saumon_Chabrier()
//                2. In_*[] and Out[] must NOT overlap
//
// Parameter   :  Mode       : To support multiple modes in this general converter
//                Out        : Output array
//                In_*       : Input array
//                AuxArray_* : Auxiliary arrays (see the Note above)
//                Table      : EoS tables
//
// Return      :  Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void EoS_General_Saumon_Chabrier( const int Mode, real Out[], const real In_Flt[], const int In_Int[],
                                       const double AuxArray_Flt[], const int AuxArray_Int[],
                                       const real *const Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
   if ( Out          == NULL )   printf( "ERROR : Out == NULL in %s !!\n", __FUNCTION__ );
   if ( In_Flt       == NULL )   printf( "ERROR : In_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Flt == NULL )   printf( "ERROR : AuxArray_Flt == NULL in %s !!\n", __FUNCTION__ );
   if ( AuxArray_Int == NULL )   printf( "ERROR : AuxArray_Int == NULL in %s !!\n", __FUNCTION__ );
#  endif // GAMER_DEBUG


   real Dens_CGS = In_Flt[0];
   real Temp_CGS = In_Flt[1];
   real Log_Temp_CGS = LOG10( Temp_CGS );
   real cooling = 0.0;
   real heating = 0.0;
   real CoolingRate = 0.0;
   
   if ( Temp_CGS < 10035.0 )
   {
      real Electron_Density= 2.4e-3 * POW( Temp_CGS/100, 0.25) / 0.5;
      real x =  MAX( MIN( Electron_Density/Dens_CGS , 0.1 ), 3.5e-4 * 0.4 );
      real G0 = 1/1.7;
      real param = G0 * SQRT(Temp_CGS)/Dens_CGS/x;
      real epsilon = 4.9e-2 / ( 1.0 + POW( param / 1925.0, 0.73 ) ) + 3.7e-2 * POW(Temp_CGS/1e4, 0.7) / ( 1.0 + param / 5e3 );

      real Carbon_Cooling = 92 * 1.38e-16 * 2 * ( 2.8e-7 / SQRT( Temp_CGS/100 ) * x + 8e-10 * POW( Temp_CGS/100, 0.07) ) * 3.5e-4 * 0.4 * EXP(-92/Temp_CGS);
      real Oxygen_Cooling = 1e-26 * SQRT(Temp_CGS) * ( 24 * EXP(-228/Temp_CGS) + 7 * EXP(-326/Temp_CGS) ) * 4.5e-4;
      real Hydrogen_Cooling = 7.3e-19 * x * EXP(-118400/Temp_CGS);
      real MetaStable_Carbon_Cooling = 6.2e4 * 1.38e-16 * ( 2.38e-8 / SQRT(Temp_CGS/1e4) * x + 1e-12 ) * EXP(-6.2e4/Temp_CGS) * 3.5e-4 * 0.4;
      real Recombination_Cooling = 4.65e-30 * POW(Temp_CGS, 0.94 ) * POW( param, 0.74/POW(Temp_CGS, 0.068) ) * x;
      real MetaStable_Oxygen_Cooling = 0.0;
      if ( Temp_CGS <= 1e4 )
      {
         MetaStable_Oxygen_Cooling  = 2.3e4 * 1.38e-16 / 3 * ( 5.1e-9 * POW(Temp_CGS/1e4, 0.57) * x + 1e-12 ) * EXP(-2.3e4/Temp_CGS) * 4.5e-4;
         MetaStable_Oxygen_Cooling += 2.9e4 * 1.38e-16 / 3 * ( 2.5e-9 * POW(Temp_CGS/1e4, 0.57) * x + 1e-12 ) * EXP(-4.9e4/Temp_CGS) * 4.5e-4;
         MetaStable_Oxygen_Cooling += 2.6e4 * 1.38e-16 / 1 * ( 5.2e-9 * POW(Temp_CGS/1e4, 0.57) * x + 1e-12 ) * EXP(-2.6e4/Temp_CGS) * 4.5e-4;
      }
      else
      {
         MetaStable_Oxygen_Cooling  = 2.3e4 * 1.38e-16 / 3 * ( 5.1e-9 * POW(Temp_CGS/1e4, 0.17) * x + 1e-12 ) * EXP(-2.3e4/Temp_CGS) * 4.5e-4;
         MetaStable_Oxygen_Cooling += 2.9e4 * 1.38e-16 / 3 * ( 2.5e-9 * POW(Temp_CGS/1e4, 0.13) * x + 1e-12 ) * EXP(-4.9e4/Temp_CGS) * 4.5e-4;
         MetaStable_Oxygen_Cooling += 2.6e4 * 1.38e-16 / 1 * ( 5.2e-9 * POW(Temp_CGS/1e4, 0.15) * x + 1e-12 ) * EXP(-2.6e4/Temp_CGS) * 4.5e-4;
      }
      cooling = Carbon_Cooling + Oxygen_Cooling + Hydrogen_Cooling + MetaStable_Carbon_Cooling + MetaStable_Oxygen_Cooling + Recombination_Cooling;
      heating = 1e-24 * epsilon * G0;
      CoolingRate = heating * Dens_CGS - cooling * SQR(Dens_CGS);
      
   }
   else
   {
      if ( Log_Temp_CGS < 4.0 )
      {
         cooling = 0.1343 * CUBE(Log_Temp_CGS) - 1.3906 * SQR(Log_Temp_CGS) + 5.1554 * Log_Temp_CGS - 31.967;
      }
      else if ( Log_Temp_CGS < 4.25 )
      {
         cooling = 12.64 * Log_Temp_CGS - 75.56;
      }
      else if ( Log_Temp_CGS < 4.35 )
      {
         cooling = -0.3 * Log_Temp_CGS - 20.565;
      }
      else if ( Log_Temp_CGS < 4.90 )
      {
         cooling = 1.745 * Log_Temp_CGS - 29.463;
      }
      else if ( Log_Temp_CGS < 5.40 )
      {
         cooling = -20.9125;
      }
      else if ( Log_Temp_CGS < 5.90 )
      {
         cooling = -1.795 * Log_Temp_CGS - 11.219;
      }
      else if ( Log_Temp_CGS < 6.20 )
      {
         cooling = -21.8095;
      }
      else if ( Log_Temp_CGS < 6.70 )
      {
         cooling = -1.261 * Log_Temp_CGS - 13.991;
      }
      else
      {
         cooling = -22.44;
      }
      
      cooling = -1.0 * POW( 10.0, cooling );
      CoolingRate = heating * Dens_CGS + cooling * SQR(Dens_CGS);
   }

   Out[0] = CoolingRate;


} // FUNCTION : EoS_General_Saumon_Chabrier



// =============================================
// III. Set EoS initialization functions
// =============================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE EoS_DE2P_t EoS_DensEint2Pres_Ptr = EoS_DensEint2Pres_Saumon_Chabrier;
FUNC_SPACE EoS_DP2E_t EoS_DensPres2Eint_Ptr = EoS_DensPres2Eint_Saumon_Chabrier;
FUNC_SPACE EoS_DP2C_t EoS_DensPres2CSqr_Ptr = EoS_DensPres2CSqr_Saumon_Chabrier;
FUNC_SPACE EoS_DE2T_t EoS_DensEint2Temp_Ptr = EoS_DensEint2Temp_Saumon_Chabrier;
FUNC_SPACE EoS_DT2P_t EoS_DensTemp2Pres_Ptr = EoS_DensTemp2Pres_Saumon_Chabrier;
FUNC_SPACE EoS_DE2S_t EoS_DensEint2Entr_Ptr = EoS_DensEint2Entr_Saumon_Chabrier;
FUNC_SPACE EoS_GENE_t EoS_General_Ptr       = EoS_General_Saumon_Chabrier;

//-----------------------------------------------------------------------------------------
// Function    :  EoS_SetCPU/GPUFunc_Saumon_Chabrier
// Description :  Return the function pointers of the CPU/GPU EoS routines
//
// Note        :  1. Invoked by EoS_Init_Saumon_Chabrier()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      EoS_SetFunc_Saumon_Chabrier( CPU_FuncPtr, GPU_FuncPtr );
//
//                3. Call-by-reference
//
// Parameter   :  EoS_DensEint2Pres_CPU/GPUPtr : CPU/GPU function pointers to be set
//                EoS_DensPres2Eint_CPU/GPUPtr : ...
//                EoS_DensPres2CSqr_CPU/GPUPtr : ...
//                EoS_DensEint2Temp_CPU/GPUPtr : ...
//                EoS_DensTemp2Pres_CPU/GPUPtr : ...
//                EoS_DensEint2Entr_CPU/GPUPtr : ...
//                EoS_General_CPU/GPUPtr       : ...
//
// Return      :  EoS_DensEint2Pres_CPU/GPUPtr, EoS_DensPres2Eint_CPU/GPUPtr,
//                EoS_DensPres2CSqr_CPU/GPUPtr, EoS_DensEint2Temp_CPU/GPUPtr,
//                EoS_DensTemp2Pres_CPU/GPUPtr, EoS_DensEint2Entr_CPU/GPUPtr,
//                EoS_General_CPU/GPUPtr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void EoS_SetGPUFunc_Saumon_Chabrier( EoS_DE2P_t &EoS_DensEint2Pres_GPUPtr,
                                   EoS_DP2E_t &EoS_DensPres2Eint_GPUPtr,
                                   EoS_DP2C_t &EoS_DensPres2CSqr_GPUPtr,
                                   EoS_DE2T_t &EoS_DensEint2Temp_GPUPtr,
                                   EoS_DT2P_t &EoS_DensTemp2Pres_GPUPtr,
                                   EoS_DE2S_t &EoS_DensEint2Entr_GPUPtr,
                                   EoS_GENE_t &EoS_General_GPUPtr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Pres_GPUPtr, EoS_DensEint2Pres_Ptr, sizeof(EoS_DE2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2Eint_GPUPtr, EoS_DensPres2Eint_Ptr, sizeof(EoS_DP2E_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensPres2CSqr_GPUPtr, EoS_DensPres2CSqr_Ptr, sizeof(EoS_DP2C_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Temp_GPUPtr, EoS_DensEint2Temp_Ptr, sizeof(EoS_DE2T_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensTemp2Pres_GPUPtr, EoS_DensTemp2Pres_Ptr, sizeof(EoS_DT2P_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_DensEint2Entr_GPUPtr, EoS_DensEint2Entr_Ptr, sizeof(EoS_DE2S_t) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &EoS_General_GPUPtr,       EoS_General_Ptr,       sizeof(EoS_GENE_t) )  );
}

#else // #ifdef __CUDACC__

void EoS_SetCPUFunc_Saumon_Chabrier( EoS_DE2P_t &EoS_DensEint2Pres_CPUPtr,
                                   EoS_DP2E_t &EoS_DensPres2Eint_CPUPtr,
                                   EoS_DP2C_t &EoS_DensPres2CSqr_CPUPtr,
                                   EoS_DE2T_t &EoS_DensEint2Temp_CPUPtr,
                                   EoS_DT2P_t &EoS_DensTemp2Pres_CPUPtr,
                                   EoS_DE2S_t &EoS_DensEint2Entr_CPUPtr,
                                   EoS_GENE_t &EoS_General_CPUPtr )
{
   EoS_DensEint2Pres_CPUPtr = EoS_DensEint2Pres_Ptr;
   EoS_DensPres2Eint_CPUPtr = EoS_DensPres2Eint_Ptr;
   EoS_DensPres2CSqr_CPUPtr = EoS_DensPres2CSqr_Ptr;
   EoS_DensEint2Temp_CPUPtr = EoS_DensEint2Temp_Ptr;
   EoS_DensTemp2Pres_CPUPtr = EoS_DensTemp2Pres_Ptr;
   EoS_DensEint2Entr_CPUPtr = EoS_DensEint2Entr_Ptr;
   EoS_General_CPUPtr       = EoS_General_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void EoS_SetAuxArray_Saumon_Chabrier( double [], int [] );
void EoS_SetCPUFunc_Saumon_Chabrier( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_DE2S_t &, EoS_GENE_t & );
#ifdef GPU
void EoS_SetGPUFunc_Saumon_Chabrier( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t &, EoS_DE2T_t &, EoS_DT2P_t &, EoS_DE2S_t &, EoS_GENE_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  EoS_Init_Saumon_Chabrier
// Description :  Initialize EoS
//
// Note        :  1. Set auxiliary arrays by invoking EoS_SetAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU EoS routines by invoking EoS_SetCPU/GPUFunc_*()
//                3. Invoked by EoS_Init()
//                   --> Enable it by linking to the function pointer "EoS_Init_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void EoS_Init_Saumon_Chabrier()
{
   Read_Saumon_Chabrier_Table( ISM_Saumon_Chabrier_Table );
   Fill_Hole( h_EoS_Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP], SAUMON_CHABRIER_TAB_DENSEINT2TEMP );
   Fill_Hole( h_EoS_Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES], SAUMON_CHABRIER_TAB_DENSEINT2PRES );
   Fill_Hole( h_EoS_Table[SAUMON_CHABRIER_TAB_DENSEINT2CS], SAUMON_CHABRIER_TAB_DENSEINT2CS );
   
   EoS_SetAuxArray_Saumon_Chabrier( EoS_AuxArray_Flt, EoS_AuxArray_Int );
   EoS_SetCPUFunc_Saumon_Chabrier( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr,
                                 EoS_DensPres2CSqr_CPUPtr, EoS_DensEint2Temp_CPUPtr,
                                 EoS_DensTemp2Pres_CPUPtr, EoS_DensEint2Entr_CPUPtr,
                                 EoS_General_CPUPtr );
#  ifdef GPU
   EoS_SetGPUFunc_Saumon_Chabrier( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr,
                                 EoS_DensPres2CSqr_GPUPtr, EoS_DensEint2Temp_GPUPtr,
                                 EoS_DensTemp2Pres_GPUPtr, EoS_DensEint2Entr_GPUPtr,
                                 EoS_General_GPUPtr );
#  endif

   Create_Table( SAUMON_CHABRIER_TAB_DENSTEMP2EINT ); // Dens, Temp -> Eint
   Create_Table( SAUMON_CHABRIER_TAB_DENSPRES2EINT ); // Dens, Pres -> Eint

} // FUNCTION : EoS_Init_Saumon_Chabrier

#endif // #ifndef __CUDACC__



#endif // #if ( MODEL == HYDRO && EOS == EOS_SAUMON_CHABRIER )
