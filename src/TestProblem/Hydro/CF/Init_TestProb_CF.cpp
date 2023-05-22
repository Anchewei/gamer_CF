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
double            rho_AD;                         // adiabatic density thresheld

static double     CF_n0;
static double     CF_vflow;
static double     CF_Mach;
static double     CF_B;
static double     CF_theta_B;
static char       CF_Tur_Table[MAX_STRING];
// =======================================================================================

#ifdef FEEDBACK
void FB_Init_SinkAccretion();
#endif
#  if ( EOS == EOS_USER )
void EoS_Init_Barotropic();
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



// errors
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
//                3. Must call EoS_Init() before calling any other EoS routine
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

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "CF_n0",             &CF_n0,                 1.0,           1.0,              NoMax_double      );
   ReadPara->Add( "CF_vflow",          &CF_vflow,              0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "CF_Mach",           &CF_Mach,               0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "CF_B",              &CF_B,                  0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "CF_theta_B",        &CF_theta_B,            0.0,           0.0,              NoMax_double      );
   ReadPara->Add( "CF_Tur_Table",      CF_Tur_Table,           NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "rho_AD",            &rho_AD,                1e-14,         0.0,              NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
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

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   Cs = SQRT( ( Const_kB*ISO_TEMP/UNIT_E ) / ( MOLECULAR_WEIGHT*Const_amu/UNIT_M ));
   rho_AD /= UNIT_D;
   CF_theta_B = CF_theta_B*M_PI/180; // degree to radian

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

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
      Aux_Message( stdout, "  test problem ID       = %d\n",            TESTPROB_ID                          );
      Aux_Message( stdout, "  Initial density       = %13.7e /cm3\n",   CF_n0                                );
      Aux_Message( stdout, "  Flow velocity         = %13.7e km/s\n",   CF_vflow                             );
      Aux_Message( stdout, "  Mach number           = %13.7e \n",       CF_Mach                              );
      Aux_Message( stdout, "  Temperature           = %13.7e K\n",      ISO_TEMP                             );
      Aux_Message( stdout, "  Sound speed           = %13.7e km/s\n",   Cs*UNIT_V/Const_km                   );
      Aux_Message( stdout, "  Magnetic field        = %13.7e G\n",      CF_B                                 );
      Aux_Message( stdout, "  Magnetic field angle  = %13.7e radian\n", CF_theta_B                           );
      Aux_Message( stdout, "  Turbulence table      = %s\n",            CF_Tur_Table                         );
      Aux_Message( stdout, "  rho_AD                = %13.7e g/cm3\n",  rho_AD                               );
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
void Load_Turbulence()
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
   Vrms_Scale = CF_Mach * Cs / Vrms;
} // Function : Load_Turbulence

//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
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

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
   double VelX, VelY, VelZ, dir;

   int i = (int) ( ( x / BoxSize[0] ) * size );    // turbulence box index (cude)
   int j = (int) ( ( y / BoxSize[0] ) * size );
   int k = (int) ( ( z / BoxSize[0] ) * size );
   int mk = FMOD(k, size);                          // modify the k index to repeat the cude in z direction
   int index = i * SQR(size) + j * size + mk;
   if ( i < 0 || i > size   ) Aux_Error( ERROR_INFO, "index is out of bound\n,  i = %d", i  );
   if ( j < 0 || j > size   ) Aux_Error( ERROR_INFO, "index is out of bound\n,  j = %d", j  );
   if ( mk < 0 || mk > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, mk = %d", mk );
   if ( index < 0 || index > tur_table_NBin ) Aux_Error( ERROR_INFO, "index is out of bound\n, index = %d", index );
   VelX = Vrms_Scale * ( Table_VelX[ index ] - Total_VelX / Total_Vrms_Count );
   VelY = Vrms_Scale * ( Table_VelY[ index ] - Total_VelY / Total_Vrms_Count );
   if (z > BoxCenter[2]) dir = -1.0;
   else dir = 1.0;
   VelZ = dir*(CF_vflow*Const_km/UNIT_V) + Vrms_Scale * ( Table_VelZ[ index ] - Total_VelZ / Total_Vrms_Count );

   Dens = CF_n0*MOLECULAR_WEIGHT*Const_amu/UNIT_D;
   MomX = Dens * VelX;
   MomY = Dens * VelY;
   MomZ = Dens * VelZ;
   // Eint = Dens * SQR(Cs) / ( GAMMA - 1.0 );
   Pres = EoS_DensTemp2Pres_CPUPtr( Dens, ISO_TEMP, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{
   double MagX = 0.0, MagY = 0.0, MagZ = 0.0;
   MagZ = CF_B*COS(CF_theta_B)/UNIT_B;
   MagY = CF_B*SIN(CF_theta_B)/UNIT_B;

   magnetic[MAGX] = MagX;
   magnetic[MAGY] = MagY;
   magnetic[MAGZ] = MagZ;

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  User-specified function to initialize particle attributes
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_LB_Init_RedistributeByRectangular()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                ParType       : Particle type     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data

//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction( const long NPar_ThisRank, const long NPar_AllRank,
                          real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                          real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                          real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{
   // we only use star particles, so keep here empty.
}  // FUNCTION : Par_Init_ByFunction


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CF
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CF()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   if ( OPT__INIT != INIT_BY_RESTART ) Load_Turbulence();
   Init_Function_User_Ptr            = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr     = SetBFieldIC;
#  endif
// comment out Init_ByFile_User_Ptr to use the default
// Init_ByFile_User_Ptr              = NULL; // option: OPT__INIT=3;             example: Init/Init_ByFile.cpp -> Init_ByFile_Default()
   Init_Field_User_Ptr               = NULL; // set NCOMP_PASSIVE_USER;          example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
   Flag_User_Ptr                     = NULL; // option: OPT__FLAG_USER;          example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr          = NULL; // option: OPT__DT_USER;            example: Miscellaneous/Mis_GetTimeStep_User.cpp
   Mis_UserWorkBeforeNextLevel_Ptr   = NULL; //                                  example: Miscellaneous/Mis_UserWorkBeforeNextLevel.cpp
   Mis_UserWorkBeforeNextSubstep_Ptr = NULL; //                                  example: Miscellaneous/Mis_UserWorkBeforeNextSubstep.cpp
   BC_User_Ptr                       = NULL; // option: OPT__BC_FLU_*=4;         example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
#  ifdef MHD
   BC_BField_User_Ptr                = NULL; // option: OPT__BC_FLU_*=4;
#  endif
   Flu_ResetByUser_Func_Ptr          = NULL; // option: OPT__RESET_FLUID;        example: Fluid/Flu_ResetByUser.cpp
   Init_DerivedField_User_Ptr        = NULL; // option: OPT__OUTPUT_USER_FIELD;  example: Fluid/Flu_DerivedField_User.cpp
   Output_User_Ptr                   = NULL; // option: OPT__OUTPUT_USER;        example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr               = NULL; // option: OPT__RECORD_USER;        example: Auxiliary/Aux_Record_User.cpp
   Init_User_Ptr                     = NULL; // option: none;                    example: none
   End_User_Ptr                      = NULL; // option: none;                    example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  ifdef GRAVITY
   Init_ExtAcc_Ptr                   = NULL; // option: OPT__EXT_ACC;            example: SelfGravity/CPU_Gravity/CPU_ExtAcc_PointMass.cpp
   End_ExtAcc_Ptr                    = NULL;
   Init_ExtPot_Ptr                   = NULL; // option: OPT__EXT_POT;            example: SelfGravity/CPU_Poisson/CPU_ExtPot_PointMass.cpp
   End_ExtPot_Ptr                    = NULL;
   Poi_AddExtraMassForGravity_Ptr    = NULL; // option: OPT__GRAVITY_EXTRA_MASS; example: none
   Poi_UserWorkBeforePoisson_Ptr     = NULL; // option: none;                    example: SelfGravity/Poi_UserWorkBeforePoisson.cpp
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr           = Par_Init_ByFunction; // option: PAR_INIT=1;              example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr       = NULL; // set PAR_NATT_USER;               example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
#  endif
#  ifdef FEEDBACK
   FB_Init_User_Ptr                  = FB_Init_SinkAccretion;
#  endif
#  if ( EOS == EOS_USER )
   EoS_Init_Ptr                      = EoS_Init_Barotropic; // option: EOS in the Makefile;     example: EoS/User_Template/CPU_EoS_User_Template.cpp
   EoS_End_Ptr                       = NULL;
#  endif
#  endif // #if ( MODEL == HYDRO )
   Src_Init_User_Ptr                 = NULL; // option: SRC_USER;                example: SourceTerms/User_Template/CPU_Src_User_Template.cpp


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CF
