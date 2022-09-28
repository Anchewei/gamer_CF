#include "GAMER.h"
#include "EoS_Barotrop1D.h"

#if ( MODEL == HYDRO && EOS == EOS_BAROTROP1D )
extern int    nRho;
extern double Min_Rho;
extern double Max_Rho;
extern double dRho; 
extern double *Table_All;
extern double *Table_Rho;
extern double *Table_Temp;

void Read_Barotrop1D_Table( char *ISM_Baratrop1D_Table_Name )
{
    if ( MPI_Rank == 0 ) Aux_Message( stdout, "%s ...\n", __FUNCTION__ );
    if ( MPI_Rank == 0 ) Aux_Message( stdout, "   Reading nuclear EoS table: %s\n", ISM_Baratrop1D_Table_Name );

// check file existence
    if ( !Aux_CheckFileExist(ISM_Baratrop1D_Table_Name) ) Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", ISM_Baratrop1D_Table_Name );

    const bool RowMajor_No  = false;           // load data into the column major
    const bool AllocMem_Yes = true;            // allocate memory for ISM_Velocity_Perturbation
    const int  size = 599;
    int Table_NCol = 2;
    int Table_TargetCols[2] = {0, 1};
    int Table_ColIdx_Rho = 0;
    int Table_ColIdx_Temp = 1;
    Min_Rho = -24.000000000000000;
    Max_Rho = 6.0100200400801604;
    dRho = 5.01002004008016047e-2;

    nRho = Aux_LoadTable( Table_All, ISM_Baratrop1D_Table_Name, Table_NCol, Table_TargetCols, RowMajor_No, AllocMem_Yes );
    if ( Table_All == NULL ) Aux_Error( ERROR_INFO, "Load %s fail!\n", ISM_Baratrop1D_Table_Name );
    if ( nRho != size ) Aux_Error( ERROR_INFO, "Load result wrong, nRho = %d!\n", nRho );
    
    Table_Rho     = Table_All + Table_ColIdx_Rho * nRho;
    Table_Temp    = Table_All + Table_ColIdx_Temp * nRho;

    h_EoS_Table[BAROTROP1D_TAB_ALL] = Table_All;
    h_EoS_Table[BAROTROP1D_TAB_RHO] = Table_Rho;
    h_EoS_Table[BAROTROP1D_TAB_TEMP] = Table_Temp;

    if ( MPI_Rank == 0 ) Aux_Message( stdout, "%s ...done\n", __FUNCTION__ );

}

#endif // #if ( MODEL == HYDRO && EOS == EOS_BAROTROP1D )