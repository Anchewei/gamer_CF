#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double     Cs;                             // sound spped
static double     R0;
static double     Rho0;
static double     Omega0;
static double     Core_Mass;
static double     Delta_Dens;
static double     Dens_Contrast;
double            rho_AD_BB;                      // adiabatic density thresheld
// =======================================================================================

#ifdef FEEDBACK
void FB_Init_SinkAccretion();
#endif
#  if ( EOS == EOS_USER )
void EoS_Init_Barotropic_BBTest();
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

#  ifndef STAR_FORMATION
   Aux_Error( ERROR_INFO, "STAR_FORMATION must be enabled !!\n" );
#  endif

#  ifndef FEEDBACK
   Aux_Error( ERROR_INFO, "FEEDBACK must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif

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
   ReadPara->Add( "R0",                &R0,                    0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Omega0",            &Omega0,                0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Core_Mass",         &Core_Mass,             0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Delta_Dens",        &Delta_Dens,            0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Dens_Contrast",     &Dens_Contrast,         0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "rho_AD_BB",         &rho_AD_BB,             0.0,           0.0,              NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) set the problem-specific derived parameters

   Core_Mass *= Const_Msun;
   Cs = SQRT( ( Const_kB*ISO_TEMP/UNIT_E ) / ( MOLECULAR_WEIGHT*Const_amu/UNIT_M ));
   R0 /= UNIT_L;
   Rho0 = 3.0 * Core_Mass / (4.0 * M_PI * CUBE(R0));
   Omega0 /= 1/UNIT_T;
   rho_AD_BB /= UNIT_D;

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
      Aux_Message( stdout, "  test problem ID       = %d\n",             TESTPROB_ID                         );
      Aux_Message( stdout, "  Sound speed           = %13.7e km/s\n",   Cs*UNIT_V/Const_km                   );
      Aux_Message( stdout, "  R0                    = %13.7e cm\n",     R0*UNIT_L                            );
      Aux_Message( stdout, "  Rho0                  = %13.7e g/cm3\n",  Rho0*UNIT_D                          );
      Aux_Message( stdout, "  Omega0                = %13.7e /s\n",     Omega0*UNIT_T                        );
      Aux_Message( stdout, "  Core Mass             = %13.7e M_sun\n",  Core_Mass/Const_Msun                 );
      Aux_Message( stdout, "  rho_AD_BB             = %13.7e \n",       rho_AD_BB                            );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter

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
   const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   const double dx = x - BoxCenter[0];
   const double dy = y - BoxCenter[0];
   const double dz = z - BoxCenter[0];
   const double Rs = SQRT( SQR(dx) + SQR(dy) + SQR(dz) );

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
   double VelX, VelY, VelZ;

   if ( Rs < R0 )
   {
      Dens = Rho0 * (1 + Delta_Dens * COS(2 * ATAN(dy/dx)));
      VelX -= Omega0 * dy;
      VelY += Omega0 * dx;
   }
   else
   {
      Dens = Rho0 / Dens_Contrast;
   }

   Pres = EoS_DensTemp2Pres_CPUPtr( Dens, ISO_TEMP, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );

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
void Par_Init_ByFunction_BBTest( const long NPar_ThisRank, const long NPar_AllRank,
                          real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                          real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                          real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{
   // we only use star particles, so keep here empty.
}  // FUNCTION : Par_Init_ByFunction
#  endif



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_BBTest
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_BBTest()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();

// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr = SetGridIC;
#  endif // #if ( MODEL == HYDRO )
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr           = Par_Init_ByFunction_BBTest; // option: PAR_INIT=1;              example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr       = NULL; // set PAR_NATT_USER;               example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
#  endif
#  ifdef FEEDBACK
   FB_Init_User_Ptr                  = FB_Init_SinkAccretion;
#  endif
#  if ( EOS == EOS_USER )
   EoS_Init_Ptr                      = EoS_Init_Barotropic_BBTest; // option: EOS in the Makefile;     example: EoS/User_Template/CPU_EoS_User_Template.cpp
   EoS_End_Ptr                       = NULL;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_BBTest