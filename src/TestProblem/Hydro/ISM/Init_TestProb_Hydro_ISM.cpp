#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double    *ISM_Prof = NULL;                // Radial profile of initial condition
static int        ISM_Prof_NBin;                  // Number of radial bins in the input profile
static int        ISM_NCol;                       // Number of columns read from the input profile
static int       *ISM_TargetCols = new int [6];   // Index of columns read from the input profile
static int        IMS_ColIdx_X;                   // Column index of x coordinate in the input profile
static int        IMS_ColIdx_Y;                   // Column index of y coordinate in the input profile
static int        IMS_ColIdx_Z;                   // Column index of z coordinate in the input profile
static int        ISM_ColIdx_VelX;                // Column index of x direction velocity in the input profile
static int        ISM_ColIdx_VelY;                // Column index of y direction velocity in the input profile
static int        ISM_ColIdx_VelZ;                // Column index of z direction velocity in the input profile
static double     *Table_X;
static double     *Table_Y;
static double     *Table_Z;
static double     *Table_VelX;
static double     *Table_VelY;
static double     *Table_VelZ;
static int        size;
static double     Total_VelX;
static double     Total_VelY;
static double     Total_VelZ;
static double     Total_VelX_SQR;
static double     Total_VelY_SQR;
static double     Total_VelZ_SQR;
static double     Vrms;
static double     Vrms_Scale;
static int        Total_Vrms_Count;
static double     Total_Mass;
static double     Mag_Rot_Angle_Rad;
static double     Bg_Dens_Contrast;
static double     Cs;
static double     Rho0;
static double     R0;
static double     Zeta;
static double     Omega0;
static double     Mass_Rad[1000];
static double     Res_int;
static double     Egrav2;
static double     B0;
static double     Mass_Sph;
static double     Min_Column_Dens;
static double     Max_Column_Dens;
// =======================================================================================




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
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifndef ISM
   Aux_Error( ERROR_INFO, "ISM must be enabled !!\n" );
#  endif

   if ( !OPT__UNIT ) Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#ifdef ISM
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
   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   ISM_NCol = 6;
   ISM_TargetCols[0] =  0;
   ISM_TargetCols[1] =  1;
   ISM_TargetCols[2] =  2;
   ISM_TargetCols[3] =  3;
   ISM_TargetCols[4] =  4;
   ISM_TargetCols[5] =  5;
   IMS_ColIdx_X      =  0;
   IMS_ColIdx_Y      =  1;
   IMS_ColIdx_Z      =  2;
   ISM_ColIdx_VelX   =  3;
   ISM_ColIdx_VelY   =  4;
   ISM_ColIdx_VelZ   =  5;
   Total_VelX = 0.0;
   Total_VelY = 0.0;
   Total_VelZ = 0.0;
   Total_VelX_SQR = 0.0;
   Total_VelY_SQR = 0.0;
   Total_VelZ_SQR = 0.0;
   Vrms = 0.0;
   Vrms_Scale = 0.0;
   Total_Vrms_Count = 0;
   size = 100;
   Egrav2 = 0.0;
   Bg_Dens_Contrast = 10.0;


// (1-3) check the runtime parameters
#  if ( EOS != EOS_ISOTHERMAL && EOS != EOS_SAUMON_CHABRIER && EOS != EOS_BAROTROP1D && EOS != EOS_ANALYTICAL_BAROTROP)
   Aux_Error( ERROR_INFO, "Must enable EOS_ISOTHERMAL or EOS_Saumon_Chabrier or EOS_Barotrop1D or EOS_ANALYTICAL_BAROTROP\n" );
#  endif


// (2) set the problem-specific derived parameters
   ISM_UNIT_T2 = SQR(MOLECULAR_WEIGHT) * SQR(Const_amu) * SQR(Const_pc) * Const_NewtonG / Const_kB;
   switch (ISM_Prob)
   {
   case 0: // Boss_Bodenheimer 
      Mag_Rot_Angle_Rad = ISM_Mag_Rot_Angle * 180.0 / M_PI;
      Cs = SQRT( ISM_Bg_Temp / ISM_UNIT_T2 );
      R0 = ISM_Alpha * 2 * Const_NewtonG * ISM_Core_Mass * MOLECULAR_WEIGHT * Const_amu / (5 * Const_kB * ISM_Bg_Temp) * UNIT_M / UNIT_L;
      Rho0 = 3.0 * ISM_Core_Mass / (4.0 * M_PI * CUBE(R0));
      Omega0 = SQRT( ISM_Beta * 4.0 * M_PI * Rho0 );
      B0 = SQRT(4.0 * M_PI /5.0)/0.53*(ISM_Crit*Rho0*R0);
      break;

   case 1:  // Bonnor Ebert like
      Mag_Rot_Angle_Rad = ISM_Mag_Rot_Angle * 180.0 / M_PI;
      Zeta = SQRT( ISM_Dens_Contrast - 1.0 );

      for ( int i = 0; i < 1000; i++ )
      {
         Res_int += LOG( 1.0 + SQR(Zeta / 1000 * (i + 1) ) ) * Zeta / 1000;
         Mass_Rad[i] = (i + 1) * Zeta / 1000 * LOG( 1.0 + SQR(Zeta / 1000 * (i + 1) ) ) - Res_int;
         Egrav2 += ( (i + 1) * Zeta / 1000)/( 1.0 + SQR(Zeta/1000*(i + 1))) * Zeta/1000 * Mass_Rad[i];
      }
      Res_int = Zeta * LOG( 1.0 + SQR(Zeta)) - Res_int;

      Cs = SQRT( ISM_Bg_Temp / ISM_UNIT_T2 );
      R0 = ISM_Core_Mass / (2.0 * M_PI * ISM_Axis_Ratio * Res_int) * SQR(ISM_FF_SCT) / (3.0 * M_PI / 32.0) / SQR(Cs);
      Rho0 = ISM_Core_Mass / (2.0 * M_PI * ISM_Axis_Ratio * Res_int) / CUBE(R0);
      Omega0 = ISM_FF_RT * 2.0 * M_PI * SQRT( 32.0 * Rho0 / 3.0 / M_PI );
      Egrav2 *= 8.0 * SQR( M_PI )* SQR(Rho0)* POW( R0, 5);

      if ( ISM_Uniform_BField ) B0 = ISM_Core_Mass * ISM_Crit/( M_PI * SQR( R0 * Zeta) ) / ( SQRT(5.0) /(3.0 *M_PI)*0.53)/SQRT(4.0*M_PI);
      else B0 = ISM_FF_ACT * SQRT(32.0/3.0/M_PI)*Rho0*R0;

      Mass_Sph = Rho0 / ISM_Dens_Contrast * amr->dh[3];
      Min_Column_Dens = amr->BoxSize[0] * Rho0 / ISM_Dens_Contrast / Bg_Dens_Contrast;
      Max_Column_Dens = Rho0 * Rho0 * ATAN(Zeta) + (amr->BoxSize[0] - 2.0 * R0 *Zeta) * Rho0 / ISM_Dens_Contrast / Bg_Dens_Contrast;

      break;

   default:
      Aux_Error( ERROR_INFO, "ISM_Prob = %d is not valid\n", ISM_Prob );
   }


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const real   ISM_Free_Fall_Time =  SQRT(3.0 * M_PI / 32.0 / Const_NewtonG / Rho0 / UNIT_D ) / UNIT_T;
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 3 * ISM_Free_Fall_Time;

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
      Aux_Message( stdout, "  TESTPROB_ID               = %d\n",      TESTPROB_ID                            );
      Aux_Message( stdout, "  ISM_Prob                  = %d\n",      ISM_Prob                               );
      Aux_Message( stdout, "  ISM_Core_Mass             = %20.14e\n", ISM_Core_Mass                          );
      Aux_Message( stdout, "  ISM_Axis_Ratio            = %20.14e\n", ISM_Axis_Ratio                         );
      Aux_Message( stdout, "  ISM_Dens_Contrast         = %20.14e\n", ISM_Dens_Contrast                      );
      Aux_Message( stdout, "  ISM_Alpha                 = %20.14e\n", ISM_Alpha                              );
      Aux_Message( stdout, "  ISM_Beta                  = %20.14e\n", ISM_Beta                               );
      Aux_Message( stdout, "  ISM_Delta_Dens            = %20.14e\n", ISM_Delta_Dens                         );
      Aux_Message( stdout, "  ISM_Crit                  = %20.14e\n", ISM_Crit                               );
      Aux_Message( stdout, "  ISM_Mach                  = %20.14e\n", ISM_Mach                               );
      Aux_Message( stdout, "  ISM_FF_SCT                = %20.14e\n", ISM_FF_SCT                             );
      Aux_Message( stdout, "  ISM_FF_RT                 = %20.14e\n", ISM_FF_RT                              );
      Aux_Message( stdout, "  ISM_FF_ACT                = %20.14e\n", ISM_FF_ACT                             );
      Aux_Message( stdout, "  ISM_FF_VCT                = %20.14e\n", ISM_FF_VCT                             );
      Aux_Message( stdout, "  ISM_Mag_Rot_Angle         = %20.14e\n", ISM_Mag_Rot_Angle                      );
      Aux_Message( stdout, "  ISM_Uniform_BField        = %20.14e\n", ISM_Uniform_BField                     );
      Aux_Message( stdout, "  ISM_Bg_Temp               = %20.14e\n", ISM_Bg_Temp                            );
      Aux_Message( stdout, "  ISM_N_Star                = %20.14e\n", ISM_N_Star                             );
      Aux_Message( stdout, "  ISM_UNIT_T2               = %20.14e\n", ISM_UNIT_T2                            );
      Aux_Message( stdout, "  ISM_Saumon_Chabrier_Table = %s\n",      ISM_Saumon_Chabrier_Table              );
      Aux_Message( stdout, "  ISM_Baratrop1D_Table      = %s\n",      ISM_Baratrop1D_Table                   );
      Aux_Message( stdout, "  ISM_Velocity_Turb_Table   = %s\n",      ISM_Velocity_Turb_Table                );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Boss_Bodenheimer_Velocity_Pertubation
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void Load_Boss_Bodenheimer_Velocity_Pertubation()
{
   const bool RowMajor_No  = false;           // load data into the column major
   const bool AllocMem_Yes = true;            // allocate memory for ISM_Velocity_Perturbation
   const double BoxSize[3]   = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
   const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };

   ISM_Prof_NBin = Aux_LoadTable( ISM_Prof, ISM_Velocity_Turb_Table, ISM_NCol, ISM_TargetCols, RowMajor_No, AllocMem_Yes );

   Table_X     = ISM_Prof + IMS_ColIdx_X * ISM_Prof_NBin;
   Table_Y     = ISM_Prof + IMS_ColIdx_Y * ISM_Prof_NBin;
   Table_Z     = ISM_Prof + IMS_ColIdx_Z * ISM_Prof_NBin;
   Table_VelX  = ISM_Prof + ISM_ColIdx_VelX * ISM_Prof_NBin;
   Table_VelY  = ISM_Prof + ISM_ColIdx_VelY * ISM_Prof_NBin;
   Table_VelZ  = ISM_Prof + ISM_ColIdx_VelZ * ISM_Prof_NBin;

   for ( int i = 0; i < ISM_Prof_NBin; i++ )
   {
      double Xi = BoxSize[0] * ( (Table_X[i] + 0.5) / size ) - BoxCenter[0];
      double Yi = BoxSize[1] * ( (Table_Y[i] + 0.5) / size ) - BoxCenter[1];
      double Zi = BoxSize[2] * ( (Table_Z[i] + 0.5) / size ) - BoxCenter[2];
      double Rs = SQRT( SQR(Xi) + SQR(Yi) + SQR(Zi) );

      if ( Rs < R0 )
      {
         Total_VelX += Table_VelX[i];
         Total_VelY += Table_VelY[i];
         Total_VelZ += Table_VelZ[i];

         Total_VelX_SQR += SQR(Table_VelX[i]);
         Total_VelY_SQR += SQR(Table_VelY[i]);
         Total_VelZ_SQR += SQR(Table_VelZ[i]);

         Total_Vrms_Count ++;
      }
   }

   // Vrms = SQRT( ( Vx^2 + Vy^2 + Vz^2 ) / N + ( Vx + Vy + Vz / N) ^ 2 )
   Vrms = SQRT( (Total_VelX_SQR + Total_VelY_SQR + Total_VelZ_SQR) / Total_Vrms_Count - 
                SQR( (Total_VelX + Total_VelY + Total_VelZ) / Total_Vrms_Count ) );
   Vrms_Scale = ISM_Mach * Cs / Vrms;

   if ( MPI_Rank == 0 )
   {   
      Aux_Message( stdout, "%s                                            ...\n", __FUNCTION__                       );
      Aux_Message( stdout, "=====================================================================================\n" );
      Aux_Message( stdout, "  Vrms for given seed                   = %20.14e\n", Vrms                               );
      Aux_Message( stdout, "  Correction factor for turbulent field = %20.14e\n", Vrms                               );
      Aux_Message( stdout, "  Alpha dense core                      = %20.14e\n", ISM_Alpha                          );
      Aux_Message( stdout, "  Beta dense core                       = %20.14e\n", ISM_Beta                           );
      Aux_Message( stdout, "  Mass (Msun)                           = %20.14e\n", ISM_Core_Mass * UNIT_M / Const_Msun);
      Aux_Message( stdout, "  Rho0 (g/cm^3)                         = %20.14e\n", Rho0 * UNIT_D                      );
      Aux_Message( stdout, "  Turbulent Mach                        = %20.14e\n", ISM_Mach                           );
      Aux_Message( stdout, "  R0                                    = %20.14e\n", R0                                 );
      Aux_Message( stdout, "  BoxSize                               = %20.14e\n", BoxSize[0]                         );
      Aux_Message( stdout, "=====================================================================================\n" );
      Aux_Message( stdout, "%s                                        ...done\n", __FUNCTION__                       );
   }

} // Function : Load_Boss_Bodenheimer_Velocity_Pertubation



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Bonner_Ebert_Like_Velocity_Pertubation
// Description :  Load velocity pertubation data from file and calculatet the Vrms
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void Load_Bonner_Ebert_Like_Velocity_Pertubation()
{
   const bool RowMajor_No  = false;           // load data into the column major
   const bool AllocMem_Yes = true;            // allocate memory for ISM_Velocity_Perturbation
   const double BoxSize[3]   = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
   const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   double Eturb = 0.0, Ethermal = 0.0, Egrav = 0.0, Erot = 0.0;

   ISM_Prof_NBin = Aux_LoadTable( ISM_Prof, ISM_Velocity_Turb_Table, ISM_NCol, ISM_TargetCols, RowMajor_No, AllocMem_Yes );

   if ( ISM_Prof == NULL ) Aux_Error( ERROR_INFO, "Load %s fail!\n", ISM_Velocity_Turb_Table );
   if ( ISM_Prof_NBin != CUBE(size) ) Aux_Error( ERROR_INFO, "Load result wrong, ISM_Prof_NBin = %d!\n", ISM_Prof_NBin );

   Table_X     = ISM_Prof + IMS_ColIdx_X * ISM_Prof_NBin;
   Table_Y     = ISM_Prof + IMS_ColIdx_Y * ISM_Prof_NBin;
   Table_Z     = ISM_Prof + IMS_ColIdx_Z * ISM_Prof_NBin;
   Table_VelX  = ISM_Prof + ISM_ColIdx_VelX * ISM_Prof_NBin;
   Table_VelY  = ISM_Prof + ISM_ColIdx_VelY * ISM_Prof_NBin;
   Table_VelZ  = ISM_Prof + ISM_ColIdx_VelZ * ISM_Prof_NBin;

   for ( int i = 0; i < ISM_Prof_NBin; i++ )
   {
      double Xi = BoxSize[0] * ( (Table_X[i] + 0.5) / size) - BoxCenter[0];
      double Yi = BoxSize[1] * ( (Table_Y[i] + 0.5) / size) - BoxCenter[1];
      double Zi = BoxSize[2] * ( (Table_Z[i] + 0.5) / size) - BoxCenter[2];
      double Rs = SQR(Xi/R0) + SQR(Yi/R0) + SQR(Zi/R0/ISM_Axis_Ratio) ;

      if ( Rs < SQR( Zeta ) )
      {
         Total_VelX += ( Rho0 / ( 1.0 + Rs ) * Table_VelX[i] );
         Total_VelY += ( Rho0 / ( 1.0 + Rs ) * Table_VelY[i] );
         Total_VelZ += ( Rho0 / ( 1.0 + Rs ) * Table_VelZ[i] );

         Total_VelX_SQR += ( Rho0 / ( 1.0 + Rs ) * SQR( Table_VelX[i] ) );
         Total_VelY_SQR += ( Rho0 / ( 1.0 + Rs ) * SQR( Table_VelY[i] ) );
         Total_VelZ_SQR += ( Rho0 / ( 1.0 + Rs ) * SQR( Table_VelZ[i] ) );

         Total_Mass += ( Rho0 / (1.0 + Rs) );
         Eturb += ( Rho0/ ( 1.0 + Rs ) * ( SQR( Table_VelX[i] ) + SQR( Table_VelY[i] ) + SQR( Table_VelZ[i] ) ) );
         Erot += ( Rho0 / ( 1.0 + Rs ) * SQR(Omega0)*( SQR(Yi) +SQR(Zi) ) );
      }
   }

   Ethermal = 3.0 / 2.0 * ISM_Core_Mass * SQR( Cs );
   Egrav = 3.0/5.0 * SQR(ISM_Core_Mass) / ( R0 * Zeta );
   Erot *= 0.5 * CUBE( BoxSize[0] / size );
   Eturb *= 0.5 * CUBE( BoxSize[0] / size );

   Total_VelX /= Total_Mass;
   Total_VelY /= Total_Mass;
   Total_VelZ /= Total_Mass;
   Total_VelX_SQR /= Total_Mass;
   Total_VelY_SQR /= Total_Mass;
   Total_VelZ_SQR /= Total_Mass;

   Vrms = SQRT( Total_VelX_SQR - SQR(Total_VelX) + Total_VelY_SQR - SQR(Total_VelY) + Total_VelZ_SQR - SQR(Total_VelZ) );
   Vrms_Scale = ISM_FF_VCT * SQRT( 32.0 * Rho0 / 3.0 / M_PI ) * R0 / Vrms;
   Total_Mass *= CUBE( BoxSize[0] / size );

   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "%s                                       ...start\n", __FUNCTION__                       );
      Aux_Message( stdout, "=====================================================================================\n" );
      Aux_Message( stdout, "  Vrms                                  = %20.14e\n", Vrms                               );
      Aux_Message( stdout, "  Vrms (rescale)                        = %20.14e\n", Vrms_Scale                         );
      Aux_Message( stdout, "  ISM Core Mass                         = %20.14e\n", ISM_Core_Mass                      );
      Aux_Message( stdout, "  Total Mass                            = %20.14e\n", Total_Mass                         );
      Aux_Message( stdout, "  Central Density                       = %20.14e\n", Rho0                               );
      Aux_Message( stdout, "  Inner Radius (pc)                     = %20.14e\n", R0                                 );
      Aux_Message( stdout, "  Angular Velocity (Hz)                 = %20.14e\n", Omega0 / ( 2.0* M_PI ) / UNIT_T    );
      Aux_Message( stdout, "  Thermal Energy                        = %20.14e\n", Ethermal                           );
      Aux_Message( stdout, "  Rotational Energy                     = %20.14e\n", Erot                               );
      Aux_Message( stdout, "  Turbulent Energy                      = %20.14e\n", Eturb                              );
      Aux_Message( stdout, "  Gravitaional Energy                   = %20.14e\n", Egrav                              );
      Aux_Message( stdout, "  Gravitaional Energy (intergrate)      = %20.14e\n", Egrav2                             );
      Aux_Message( stdout, "  Thermal / gravitational energy        = %20.14e\n", Ethermal / Egrav                   );
      Aux_Message( stdout, "  Thermal / gravitational energy        = %20.14e\n", Ethermal / Egrav2                  );
      Aux_Message( stdout, "  Rotatinal / gravitational energy      = %20.14e\n", Erot / Egrav2                      );
      Aux_Message( stdout, "  Turbulent / gravitational energy      = %20.14e\n", Eturb * SQR(Vrms_Scale)/ Egrav2    );
      Aux_Message( stdout, "=====================================================================================\n" );
      Aux_Message( stdout, "%s                                         ...end\n", __FUNCTION__                       );
   }


} // Function : Load_Bonner_Ebert_Like_Velocity_Pertubation



//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Boss_Bodenheimer_GridIC
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
void Set_Boss_Bodenheimer_GridIC( real fluid[], const double x, const double y, const double z, const double Time,
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

   if ( ISM_Mach != 0.0 )
   {
      int i = (int) ( ( x / BoxSize[0] ) * size );
      int j = (int) ( ( y / BoxSize[0] ) * size );
      int k = (int) ( ( z / BoxSize[0] ) * size );
      int index = i * SQR(size) + j * size + k;
      if ( i < 0 || i > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, i = %d", i );
      if ( j < 0 || j > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, j = %d", j );
      if ( k < 0 || k > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, k = %d", k );
      if ( index < 0 || index > ISM_Prof_NBin ) Aux_Error( ERROR_INFO, "index is out of bound\n, index = %d", index );
      VelX = Vrms_Scale * ( Table_VelX[ index ] - Total_VelX / Total_Vrms_Count );
      VelY = Vrms_Scale * ( Table_VelY[ index ] - Total_VelY / Total_Vrms_Count );
      VelZ = Vrms_Scale * ( Table_VelZ[ index ] - Total_VelZ / Total_Vrms_Count );
   }

   if ( Rs < R0 )
   {
      Dens = Rho0 * (1 + ISM_Delta_Dens * COS(2 * ATAN(dy / (COS(Mag_Rot_Angle_Rad) * dx - SIN(Mag_Rot_Angle_Rad) * dz))));
      VelX += ( Omega0 * dy * COS(Mag_Rot_Angle_Rad) );
      VelY += ( Omega0 * dx * (COS(Mag_Rot_Angle_Rad) - SIN(Mag_Rot_Angle_Rad)) );
      VelZ += ( Omega0 * dy * SIN(Mag_Rot_Angle_Rad) );
   }
   else
   {
      Dens = Rho0 / ISM_Dens_Contrast;
   }

#  if ( EOS == EOS_SAUMON_CHABRIER )
   Pres = EoS_DensTemp2Pres_CPUPtr( Dens, ISM_Bg_Temp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#  else
   Eint = Dens * SQR(Cs) / ( GAMMA - 1.0 );
#  endif

   MomX = Dens * VelX;
   MomY = Dens * VelY;
   MomZ = Dens * VelZ;
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : Set_Boss_Bodenheimer_GridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Bonner_Ebert_Like_GridIC
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
void Set_Bonner_Ebert_Like_GridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double BoxSize[3]   = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
   const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
   double VelX = 0.0, VelY = 0.0, VelZ = 0.0;

   const double dx = x - BoxCenter[0];
   const double dy = y - BoxCenter[1];
   const double dz = z - BoxCenter[2];
   const double Rs = SQR(dx/R0) + SQR(dy/R0) + SQR( dz / R0 / ISM_Axis_Ratio );

   if ( Rs > SQR(Zeta) )
   {
      Dens = Rho0 / ISM_Dens_Contrast / Bg_Dens_Contrast;
#     if ( EOS == EOS_SAUMON_CHABRIER )
      Pres = EoS_DensTemp2Pres_CPUPtr( Dens, ISM_Bg_Temp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#     else
      Eint = Dens * SQR(Cs) / ( GAMMA - 1.0 );
#     endif
   }
   else
   {
      Dens = Rho0 / ( 1 + Rs );
#     if ( EOS == EOS_SAUMON_CHABRIER )
      Pres = EoS_DensTemp2Pres_CPUPtr( Dens, ISM_Bg_Temp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#     else
      Eint = Dens * SQR(Cs) / ( GAMMA - 1.0 );
#     endif
   }

   if ( Vrms_Scale != 0.0 )
   {
      int i = (int) ( ( x / BoxSize[0] ) * size );
      int j = (int) ( ( y / BoxSize[1] ) * size );
      int k = (int) ( ( z / BoxSize[2] ) * size );
      int index = i * SQR(size) + j * size + k;

      if ( i < 0 || i > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, i = %d", i );
      if ( j < 0 || j > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, j = %d", j );
      if ( k < 0 || k > size ) Aux_Error( ERROR_INFO, "index is out of bound\n, k = %d", k );
      if ( index < 0 || index > ISM_Prof_NBin ) Aux_Error( ERROR_INFO, "index is out of bound\n, index = %d", index );
   
      VelX = Vrms_Scale * ( Table_VelX[ index ] - Total_VelX );
      VelY = Vrms_Scale * ( Table_VelY[ index ] - Total_VelY );
      VelZ = Vrms_Scale * ( Table_VelZ[ index ] - Total_VelZ );
   }
   
   double R = SQR(dx/R0) + SQR(dy/R0) + SQR(dz/R0);
   if ( R < SQR(Zeta * ISM_Axis_Ratio) )
   {
      VelX -= ( Omega0 * dz * SIN(Mag_Rot_Angle_Rad) );
      VelX += ( Omega0 * dz * COS(Mag_Rot_Angle_Rad) );
      VelX += ( Omega0 * ( dx * SIN(Mag_Rot_Angle_Rad) - dy * COS(Mag_Rot_Angle_Rad) ) );
   }

   MomX = Dens * VelX;
   MomY = Dens * VelY;
   MomZ = Dens * VelZ;
   
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : Set_Bonner_Ebert_Like_GridIC

#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Boss_Bodenheimer_BFieldIC
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
void Set_Boss_Bodenheimer_BFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{
   const double BoxSize[3]   = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
   const double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
   double dx = x - BoxCenter[0];
   double dy = y - BoxCenter[1];
   double dz = z - BoxCenter[2];
   double Rc = SQRT( SQR(dx) + SQR(dy) );
   double MagX = 0.0, MagY = 0.0, MagZ = 0.0;

   if ( Rc < R0 ) MagZ = B0;
   else MagZ = B0/( POW( ISM_Dens_Contrast, 2.0/3.0) );

   magnetic[MAGX] = MagX;
   magnetic[MAGY] = MagY;
   magnetic[MAGZ] = MagZ;

} // FUNCTION : Set_Boss_Bodenheimer_BFieldIC

//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Bonner_Ebert_Like_BFieldIC
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
void Set_Bonner_Ebert_Like_BFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{
   if ( ISM_Uniform_BField )
   {
      double MagX = B0, MagY = 0.0, MagZ = 0.0;
      magnetic[MAGX] = MagX;
      magnetic[MAGY] = MagY;
      magnetic[MAGZ] = MagZ;
   }
   else
   {
      double MagX = 0.0, MagY = 0.0, MagZ = 0.0, Column_Dens = 0.0;
      double BoxSize[3] = { amr->BoxSize[0], amr->BoxSize[1], amr->BoxSize[2] };
      double BoxCenter[3] = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
      double dx = x - BoxCenter[0];
      double dy = y - BoxCenter[1];
      double dz = z - BoxCenter[2];
      double dh = amr->dh[lv];
      double Min_dh = amr->dh[3];

      if ( dh < Min_dh )
      {
         Aux_Error( ERROR_INFO, "Min_dh too large, dx = %20.14e, Min_dh = %20.14e\n", dh / amr->BoxSize[0], Min_dh / amr->BoxSize[0]);
      }

      double N = dh / Min_dh;
      double xl = dx - 0.5*dh;
      double yl = dy - 0.5*dh;
      double zl = dz - 0.5*dh;

      for ( int j = 0 ; j < N ; j++ )
      {
         for ( int k = 0; k < N; k++ )
         {
            double yy = yl + ( j + 0.5 )* Min_dh;
            double zz = zl + ( k + 0.5 )* Min_dh;
            double R = SQR(yy/R0) + SQR(zz/R0/ISM_Axis_Ratio);

            if ( R < SQR(Zeta) )
            {
               Column_Dens = R0 * Rho0 / SQRT(1.0+R) * ATAN(SQRT( (SQR(Zeta) - R)/(1.0 + R ) ) );
               Column_Dens = MAX(Column_Dens, Min_Column_Dens);
            }
            else
            {
               Column_Dens = Min_Column_Dens;
            }
            MagX += B0 * Column_Dens/Max_Column_Dens;
         }
      }

      magnetic[MAGX] = MagX / SQR(N);
      magnetic[MAGY] = MagY;
      magnetic[MAGZ] = MagZ;
   }
}

#endif // #ifdef MHD
#endif // #ifdef ISM



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ISM
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ISM()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();

#  ifdef ISM
// set the problem-specific runtime parameters
   SetParameter();

   switch (ISM_Prob)
   {
   case 0:
      if ( OPT__INIT != INIT_BY_RESTART ) Load_Boss_Bodenheimer_Velocity_Pertubation();
      Init_Function_User_Ptr         = Set_Boss_Bodenheimer_GridIC;
#     ifdef MHD
      Init_Function_BField_User_Ptr  = Set_Boss_Bodenheimer_BFieldIC;
#     endif

      break;
   
   case 1:
      if ( OPT__INIT != INIT_BY_RESTART ) Load_Bonner_Ebert_Like_Velocity_Pertubation();
      Init_Function_User_Ptr         = Set_Bonner_Ebert_Like_GridIC;
#     ifdef MHD
      Init_Function_BField_User_Ptr  = Set_Bonner_Ebert_Like_BFieldIC;
#     endif

      break;
   
   default:
      Aux_Error( ERROR_INFO, "ISM_Prob = %d is not valid\n", ISM_Prob );
   }
#  endif // ifdef ISM

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ISM
