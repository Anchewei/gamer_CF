#include "EoS_Saumon_Chabrier.h"
#include "GAMER.h"

#if (MODEL == HYDRO && EOS == EOS_SAUMON_CHABRIER)

extern int nRho;
extern int nTemp;
extern int nPres;
extern double MaxTemp;
extern double MinTemp;
extern double MaxPres;
extern double MinPres;
extern double dTemp;
extern double dPres;
extern double *Table_LogRho;
extern double *Table_DensTemp2Eint;
extern double *Table_DensPres2Eint;


void Construct_DensTemp2Eint_Table()
{
    if (MPI_Rank == 0)  Aux_Message(stdout, "%s ... \n", __FUNCTION__);

    Table_DensTemp2Eint = (double *)malloc(sizeof(double) * nRho * nTemp);
    memset(Table_DensTemp2Eint, 0, sizeof(double) * nRho * nTemp);

    for (int ir = 1; ir < nRho - 1; ir++)
    {
        for (int it = 0; it < nTemp; it++)
        {
            double Rho0 = POW(10.0, Table_LogRho[1 * nRho + ir]);
            double Temp0 = POW(10.0, (LOG10(MinTemp) + it * dTemp));
            double Eint_Old = Rho0 * Const_kB * Temp0 / (MOLECULAR_WEIGHT * Const_amu * (GAMMA - 1.0));
            if (it > 0) Eint_Old = MAX(Eint_Old, Table_DensTemp2Eint[ ( it - 1 )* nRho + ir]);

            for (int ii = 0; ii < 1000; ii++)
            {
                double Temp1 = EoS_DensEint2Temp_CPUPtr(Rho0/UNIT_D, Eint_Old/(UNIT_D*SQR(UNIT_V)), NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table);
                if (Temp1 == 0.0)
                {
                    Eint_Old = 0.0;
                    break;
                }

                double Temp2 = EoS_DensEint2Temp_CPUPtr(Rho0/UNIT_D, Eint_Old*1.001/(UNIT_D*SQR(UNIT_V)), NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table);
                if (Temp2 == 0.0)
                {
                    Eint_Old = 0.0;
                    break;
                }

                double Eint_New;
                if (FABS(Temp2 - Temp1) != 0.0)
                    Eint_New = Eint_Old - (Temp1 - Temp0) / ((Temp2 - Temp1) / (0.001 * Eint_Old));
                else
                    Eint_New = Eint_Old;

                double Epsilon = FABS(Eint_New - Eint_Old) / Eint_Old;
                Eint_Old = Eint_New;

                if (FABS(Epsilon) < 1e-4)
                    break;
                else if (ii == 999)
                {
                    if (MPI_Rank == 0) Aux_Message(stdout, "Newton for Eint(Dens,Temp) did not converge\n");
                }
            }
            Table_DensTemp2Eint[it * nRho + ir] = Eint_Old;
        }
    }
    h_EoS_Table[SAUMON_CHABRIER_TAB_DENSTEMP2EINT] = Table_DensTemp2Eint;
    if (MPI_Rank == 0)
        Aux_Message(stdout, "%s ... done\n", __FUNCTION__);
}

void Construct_DensPres2Eint_Table()
{
    if (MPI_Rank == 0) Aux_Message(stdout, "%s ... \n", __FUNCTION__);
    double Pres1, Pres1_CGS, Pres2, Pres2_CGS, Eint_New;

    Table_DensPres2Eint = (double *)malloc(sizeof(double) * nRho * nPres);
    memset(Table_DensPres2Eint, 0, sizeof(double) * nRho * nPres);

    for (int ir = 1; ir < nRho - 1; ir++)
    {
        for (int ip = 0; ip < nPres; ip++)
        {
            double Rho0 = POW(10.0, Table_LogRho[1 * nRho + ir]);
            double Pres0 = POW(10.0, ( LOG10(MinPres) + ip * dPres)) * Rho0;
            double Eint_Old = Pres0 / (GAMMA - 1.0);
            
            if (ip > 0)
			{
				if( Table_DensPres2Eint[(ip - 1) * nRho + ir] != 0.0 )
				{
					Eint_Old = Table_DensPres2Eint[(ip - 1) * nRho + ir];
				}
			}
            
            for (int ii = 0; ii < 1000; ii++)
            {   
                Pres1 = EoS_DensEint2Pres_CPUPtr(Rho0/UNIT_D, Eint_Old/(UNIT_D * SQR(UNIT_V)), NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table);
                if (Pres1 == 0.0)
                {
                    Eint_Old = 0.0;
                    break;
                }

                Pres2 = EoS_DensEint2Pres_CPUPtr(Rho0/UNIT_D, Eint_Old * 1.001/(UNIT_D * SQR(UNIT_V)), NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table);
                if (Pres2 == 0.0)
                {
                    Eint_Old = 0.0;
                    break;
                }
                
                if (FABS(Pres2 - Pres1) != 0.0)
                {
                    Pres1_CGS = Pres1 * (UNIT_D * SQR(UNIT_V));
                    Pres2_CGS = Pres2 * (UNIT_D * SQR(UNIT_V));
                    Eint_New = Eint_Old - (Pres1_CGS - Pres0 ) / ( (Pres2_CGS - Pres1_CGS ) / (0.001 * Eint_Old ));
                }
                else
                {
                    Eint_New = Eint_Old;
                }
                
                double Epsilon = FABS(Eint_New - Eint_Old) / Eint_Old;
                Eint_Old = Eint_New;

                if (FABS(Epsilon) < 1e-4)
                    break;
                else if (ii == 999)
                {
                    if (MPI_Rank == 0)
                    {
                        Aux_Message(stdout, "Newton for Eint(Dens,Pres) did not converge\n");
                    }
                }
            }

            Table_DensPres2Eint[ip * nRho + ir] = Eint_Old;
        }
    }
    h_EoS_Table[SAUMON_CHABRIER_TAB_DENSPRES2EINT] = Table_DensPres2Eint;
    if (MPI_Rank == 0)
        Aux_Message(stdout, "%s ... done\n", __FUNCTION__);
}

void Create_Table(const int Target_Table_Index )
{
    if (MPI_Rank == 0)
    {
        Aux_Message( stdout, "%s                                         ...\n", __FUNCTION__                       );
        Aux_Message( stdout, "=====================================================================================\n" );
        Aux_Message( stdout, "  NRHO                        = %d\n",  EoS_AuxArray_Int[SAUMON_CHABRIER_AUX_NRHO]       );
        Aux_Message( stdout, "  NENER                       = %d\n",  EoS_AuxArray_Int[SAUMON_CHABRIER_AUX_NENER]      );
        Aux_Message( stdout, "  NTEMP                       = %d\n",  EoS_AuxArray_Int[SAUMON_CHABRIER_AUX_NTEMP]      );
        Aux_Message( stdout, "  NPRES                       = %d\n",  EoS_AuxArray_Int[SAUMON_CHABRIER_AUX_NPRES]      );
        Aux_Message( stdout, "  MAXRHO                      = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_MAXRHO] );
        Aux_Message( stdout, "  MINRHO                      = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_MINRHO] );
        Aux_Message( stdout, "  MAXENER                     = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_MAXENER]);
        Aux_Message( stdout, "  MINENER                     = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_MINENER]);
        Aux_Message( stdout, "  MAXTEMP                     = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_MAXTEMP]);
        Aux_Message( stdout, "  MINTEMP                     = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_MINTEMP]);
        Aux_Message( stdout, "  MAXPRES                     = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_MAXPRES]);
        Aux_Message( stdout, "  MINPRES                     = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_MINPRES]);
        Aux_Message( stdout, "  DRHO                        = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_DRHO]   );
        Aux_Message( stdout, "  DENER                       = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_DENER]  );
        Aux_Message( stdout, "  DTEMP                       = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_DTEMP]  );
        Aux_Message( stdout, "  DPRES                       = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_DPRES]  );
        Aux_Message( stdout, "  UNIT_D                      = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_D] );
        Aux_Message( stdout, "  UNIT_V                      = %20.14e\n", EoS_AuxArray_Flt[SAUMON_CHABRIER_AUX_UNIT_V] );
        Aux_Message( stdout, "=====================================================================================\n" );
    }

    switch (Target_Table_Index)
    {
    case SAUMON_CHABRIER_TAB_DENSTEMP2EINT:
        Construct_DensTemp2Eint_Table();
        break;

    case SAUMON_CHABRIER_TAB_DENSPRES2EINT:
        Construct_DensPres2Eint_Table();
        break;

    default:
        Aux_Error(ERROR_INFO, "Not support construct target table = %d", Target_Table_Index);
        break;
    }

    if (MPI_Rank == 0)
        Aux_Message(stdout, "%s ... done\n", __FUNCTION__);
}

#endif // if ( MODEL == HYDRO && EOS == EOS_SAUMON_CHABRIER )
