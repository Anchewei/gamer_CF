#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double     *tur_table = NULL;              // used to store turbulence (1D)
static int        tur_table_NBin;                 // number of row in turbulence table obtained by Aux_LoadTable
static int        tur_table_Ncol;                 // number of column in turbulence table (set by user)
static int        *CF_TargetCols = new int [6];    // Index of columns read from the turbulence table 
static int        CF_ColIdx_X;                    // Column index of x coordinate in the turbulence table
static int        CF_ColIdx_Y;                    // Column index of y coordinate in the turbulence table 
static int        CF_ColIdx_Z;                    // Column index of z coordinate in the turbulence table 
static int        CF_ColIdx_VelX;                 // Column index of x direction velocity in the turbulence table 
static int        CF_ColIdx_VelY;                 // Column index of y direction velocity in the turbulence table 
static int        CF_ColIdx_VelZ;                 // Column index of z direction velocity in the turbulence table 
static double     *Table_X;                       // used to store the readed data
static double     *Table_Y;
static double     *Table_Z;
static double     *Table_VelX;
static double     *Table_VelY;
static double     *Table_VelZ;
static int        size;                           // turbulence cell number
static double     Total_VelX;                     // used to calculate Vrms
static double     Total_VelY;
static double     Total_VelZ;
static double     Total_VelX_SQR;
static double     Total_VelY_SQR;
static double     Total_VelZ_SQR;
static double     Vrms;
static double     Vrms_Scale;                     // used to rescale velocity
static int        Total_Vrms_Count;

static double     Cs;                             // sound spped

static char       CF_Tur_Table[MAX_STRING];
static double     CF_vflow;
double            rho_AD_SC;                      // adiabatic density thresheld for barotropic EOS
double            rho_PW_SC;                      // adiabatic density thresheld for piecewise EOS
static int        Plummer_like;                   // use Plummer like density profile
static int        EOS_type;                       // 0: barotropic EOS, 1: piecewise EOS

static double     Mach;
static double     Rho0;
static double     R0;
static double     Omega0;
static double     Core_Mass;
static double     Delta_Dens;
static double     Dens_Contrast;
static double     Inner_R0;                      // for Plummer_like
static double     Rho0_P;                        // for Plummer_like
// =======================================================================================

#ifdef FEEDBACK
void FB_Init_SinkAccretion();
#endif
#  if ( EOS == EOS_USER )
void EoS_Init_Barotropic_SC();
void EoS_Init_Piecewise_SC();
#  endif


//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif

// #  ifndef COMOVING
//    Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
// #  endif

// #  ifdef PARTICLE
//    Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
// #  endif

#  ifdef GRAVITY
   if ( !OPT__SELF_GRAVITY )
   Aux_Error( ERROR_INFO, "must enable OPT__SELF_GRAVITY !!\n" );

   if ( OPT__EXT_ACC )
   Aux_Error( ERROR_INFO, "must disable OPT__EXT_ACC !!\n" );

   if ( OPT__EXT_POT )
   Aux_Error( ERROR_INFO, "must disable OPT__EXT_POT !!\n" );
#  endif

   if ( MPI_Rank == 0 )
   {
#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test\n" );
#     endif

      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : non-periodic BC for fluid ??\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT != BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : non-periodic BC for gravity ??\n" );
#     endif

      if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
         Aux_Message( stderr, "WARNING : simulation domain is not cubic ??\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "CF_Tur_Table",      CF_Tur_Table,           NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "Mach",              &Mach,                  0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "CF_vflow",          &CF_vflow,              0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "R0",                &R0,                    0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Omega0",            &Omega0,                0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Core_Mass",         &Core_Mass,             0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Delta_Dens",        &Delta_Dens,            0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Dens_Contrast",     &Dens_Contrast,         0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "EOS_type",          &EOS_type,              0,             0,                NoMax_int         );
   ReadPara->Add( "rho_AD_SC",         &rho_AD_SC,             0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "rho_PW_SC",         &rho_PW_SC,             0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Plummer_like",      &Plummer_like,            0,             0,              NoMax_int         );
   ReadPara->Add( "Inner_R0",          &Inner_R0,              0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Rho0_P",            &Rho0_P,                0.0,           0.0,              NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) set the problem-specific derived parameters
   tur_table_Ncol = 6;
   CF_TargetCols[0] =  0;
   CF_TargetCols[1] =  1;
   CF_TargetCols[2] =  2;
   CF_TargetCols[3] =  3;
   CF_TargetCols[4] =  4;
   CF_TargetCols[5] =  5;
   CF_ColIdx_X      =  0;
   CF_ColIdx_Y      =  1;
   CF_ColIdx_Z      =  2;
   CF_ColIdx_VelX   =  3;
   CF_ColIdx_VelY   =  4;
   CF_ColIdx_VelZ   =  5;
   Total_VelX = 0.0;
   Total_VelY = 0.0;
   Total_VelZ = 0.0;
   Total_VelX_SQR = 0.0;
   Total_VelY_SQR = 0.0;
   Total_VelZ_SQR = 0.0;
   Vrms = 0.0;
   Vrms_Scale = 0.0;
   Total_Vrms_Count = 0;
   size = 129;

   Core_Mass *= Const_Msun;
   Cs = SQRT( ( Const_kB*ISO_TEMP/UNIT_E ) / ( MOLECULAR_WEIGHT*Const_amu/UNIT_M ));
   R0 /= UNIT_L;

   if ( Plummer_like == 1)
   {
      Rho0 = Rho0_P/UNIT_D;
   }
   else
   {
      Rho0 = 3.0 * Core_Mass / (4.0 * M_PI * CUBE(R0));
   }


   Omega0 /= 1/UNIT_T;
   rho_AD_SC /= UNIT_D;

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 5.0e-2;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID       );
      Aux_Message( stdout, "  Mach number           = %13.7e \n",       Mach                                 );
      Aux_Message( stdout, "  Flow velocity         = %13.7e km/s\n",   CF_vflow                             );
      Aux_Message( stdout, "  Sound speed           = %13.7e km/s\n",   Cs*UNIT_V/Const_km                   );
      Aux_Message( stdout, "  Turbulence table      = %s\n",            CF_Tur_Table                         );
      Aux_Message( stdout, "  Plummer_like          = %d \n",           Plummer_like                         );
      Aux_Message( stdout, "  R0                    = %13.7e cm\n",     R0*UNIT_L                            );
      Aux_Message( stdout, "  Rho0                  = %13.7e g/cm3\n",  Rho0*UNIT_D                          );
      Aux_Message( stdout, "  Omega0                = %13.7e /s\n",     Omega0*UNIT_T                        );
      Aux_Message( stdout, "  Core Mass             = %13.7e M_sun\n",  Core_Mass/Const_Msun                 );
      Aux_Message( stdout, "  EOS_type              = %d \n",           EOS_type                             );

      if ( EOS_type == 0 )
      {
         Aux_Message( stdout, "  rho_AD_SC             = %13.7e \n",       rho_AD_SC                         );
      }
      else
      {
         Aux_Message( stdout, "  rho_PW_SC             = %13.7e \n",       rho_PW_SC                         );
      }


      if ( Plummer_like == 1)
      Aux_Message( stdout, "  Inner_R0              = %13.7e \n",       Inner_R0                             );

      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter

//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Turbulence
// Description : load turbulence and calculate Vrms_Scale
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void Load_Turbulence_SC()
{
   const bool RowMajor_No  = false;           // load data into the column major
   const bool AllocMem_Yes = true;            // allocate memory for ISM_Velocity_Perturbation
   const double BoxSize[3]   = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
   const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };

   tur_table_NBin = Aux_LoadTable( tur_table, CF_Tur_Table, tur_table_Ncol, CF_TargetCols, RowMajor_No, AllocMem_Yes );

   Table_X     = tur_table + CF_ColIdx_X * tur_table_NBin;
   Table_Y     = tur_table + CF_ColIdx_Y * tur_table_NBin;
   Table_Z     = tur_table + CF_ColIdx_Z * tur_table_NBin;
   Table_VelX  = tur_table + CF_ColIdx_VelX * tur_table_NBin;
   Table_VelY  = tur_table + CF_ColIdx_VelY * tur_table_NBin;
   Table_VelZ  = tur_table + CF_ColIdx_VelZ * tur_table_NBin;

   for ( int i = 0; i < tur_table_NBin; i++ )
   {
      Total_VelX += Table_VelX[i];
      Total_VelY += Table_VelY[i];
      Total_VelZ += Table_VelZ[i];

      Total_VelX_SQR += SQR(Table_VelX[i]);
      Total_VelY_SQR += SQR(Table_VelY[i]);
      Total_VelZ_SQR += SQR(Table_VelZ[i]);

      Total_Vrms_Count ++;
   }

   // Vrms = SQRT( ( Vx^2 + Vy^2 + Vz^2 ) / N + ( Vx + Vy + Vz / N) ^ 2 )
   Vrms = SQRT( (Total_VelX_SQR + Total_VelY_SQR + Total_VelZ_SQR) / Total_Vrms_Count - 
                SQR( (Total_VelX + Total_VelY + Total_VelZ) / Total_Vrms_Count ) );
   Vrms_Scale = Mach * Cs / Vrms;
} // Function : Load_Turbulence

//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{
   const double BoxSize[3]   = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
   const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double dx = x - BoxCenter[0];
   const double dy = y - BoxCenter[0];
   const double dz = z - BoxCenter[0];
   const double Rc = SQRT( SQR(dx) + SQR(dy) );
   const double Rs = SQRT( SQR(dx) + SQR(dy) + SQR(dz) );

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
   double VelX, VelY, VelZ;

   if ( Mach != 0.0 )
   {
      int i = (int) ( ( x / BoxSize[0] ) * size );
      int j = (int) ( ( y / BoxSize[0] ) * size );
      int k = (int) ( ( z / BoxSize[0] ) * size );
      int index = i * SQR(size) + j * size + k;
      if ( i < 0 || i > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, i = %d", i );
      if ( j < 0 || j > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, j = %d", j );
      if ( k < 0 || k > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, k = %d", k );
      if ( index < 0 || index > tur_table_NBin ) Aux_Error( ERROR_INFO, "index is out of bound\n, index = %d", index );
      VelX = Vrms_Scale * ( Table_VelX[ index ] - Total_VelX / Total_Vrms_Count );
      VelY = Vrms_Scale * ( Table_VelY[ index ] - Total_VelY / Total_Vrms_Count );
      VelZ = Vrms_Scale * ( Table_VelZ[ index ] - Total_VelZ / Total_Vrms_Count );
   }

   if ( Rs < R0 )
   {
      if ( Plummer_like == 1)
      {
         Dens = Rho0 / ( SQR(Rs/Inner_R0) + 1 );
      }
      else
      {
         Dens = Rho0 * (1 + Delta_Dens * COS(2 * ATAN(dy/dx)));
      }
      VelX -= Omega0 * dy;
      VelY += Omega0 * dx;
   }
   else
   {
      if ( Plummer_like == 1)
      {
         Dens = (Rho0 / ( SQR(R0/Inner_R0) + 1 )) / Dens_Contrast;
      }
      else
      {
         Dens = Rho0 / Dens_Contrast;
      }
   }

   Pres = EoS_DensTemp2Pres_CPUPtr( Dens, ISO_TEMP, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );
   // Eint = EoS_DensPres2Eint_CPUPtr( Dens, Dens * SQR(Cs), NULL, EoS_AuxArray_Flt,
   //                                  EoS_AuxArray_Int, h_EoS_Table );
   MomX = Dens * VelX;
   MomY = Dens * VelY;
   MomZ = Dens * VelZ;
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )

#  ifdef PARTICLE
void Par_Init_ByFunction_SC( const long NPar_ThisRank, const long NPar_AllRank,
                          real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                          real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                          real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{
   // we only use star particles, so keep here empty.
}  // FUNCTION : Par_Init_ByFunction
#  endif



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_SphericalCollapse
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_SphericalCollapse()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();
   if ( OPT__INIT != INIT_BY_RESTART ) Load_Turbulence_SC();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr = SetGridIC;
#  endif // #if ( MODEL == HYDRO )
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr           = Par_Init_ByFunction_SC; // option: PAR_INIT=1;              example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr       = NULL; // set PAR_NATT_USER;               example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
#  endif
#  ifdef FEEDBACK
   FB_Init_User_Ptr        = FB_Init_SinkAccretion;
#  endif
#  if ( EOS == EOS_USER )
   if ( EOS_type == 0 )
   {
      EoS_Init_Ptr                      = EoS_Init_Barotropic_SC; // option: EOS in the Makefile;     example: EoS/User_Template/CPU_EoS_User_Template.cpp
   }
   else
   {
      EoS_Init_Ptr                      = EoS_Init_Piecewise_SC; // option: EOS in the Makefile;     example: EoS/User_Template/CPU_EoS_User_Template.cpp
   }

   EoS_End_Ptr                       = NULL;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_SphericalCollapse