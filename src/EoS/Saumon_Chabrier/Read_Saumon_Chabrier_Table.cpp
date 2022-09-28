#include "EoS_Saumon_Chabrier.h"
#include "GAMER.h"

#if ( MODEL == HYDRO && EOS == EOS_SAUMON_CHABRIER )

extern int nRho;
extern int nEner;
extern int nTemp;
extern int nPres;
extern double MaxRho;
extern double MinRho;
extern double MaxEner;
extern double MinEner;
extern double MaxTemp;
extern double MinTemp;
extern double MaxPres;
extern double MinPres;
extern double dRho;
extern double dEner;
extern double dTemp;
extern double dPres;
extern double *Table_LogRho;
extern double *Table_LogEner;
extern double *Table_DensEint2Temp;
extern double *Table_DensEint2Pres;
extern double *Table_DensEint2Cs;

void Read_Saumon_Chabrier_Table( char *ISM_Saumon_Chabrier_Table_Name )
{
    if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... \n", __FUNCTION__ );

    FILE *fp = NULL;
    fp = fopen( ISM_Saumon_Chabrier_Table_Name, "rb");
    if ( fp == NULL ) Aux_Error(ERROR_INFO, "%s is not exist.\n", ISM_Saumon_Chabrier_Table_Name);

    fseek(fp, 4, SEEK_SET); // skip blank
    fread(&nRho, sizeof(int), 1, fp);
    fread(&nEner, sizeof(int), 1, fp);
    
    size_t sizeOfTable = nRho * nEner  * sizeof(double);
    nTemp = nEner;
    nPres = nEner;
    MaxTemp = 1e5;
    MinTemp = 3.0;
    MaxPres = 1e14;
    MinPres = 1e8;

    fseek(fp, 8, SEEK_CUR); // skip blank
    fread(&MinRho, sizeof(double), 1, fp);
    fread(&MaxRho, sizeof(double), 1, fp);
    fread(&MinEner, sizeof(double), 1, fp);
    fread(&MaxEner, sizeof(double), 1, fp);
    fseek(fp, 8, SEEK_CUR); // skip yHe, because yHe is useless;

    dRho = ( MaxRho - MinRho ) / nRho;
    dEner = ( MaxEner - MinEner ) / nEner;
    dTemp = ( LOG10(MaxTemp) - LOG10(MinTemp) ) / nTemp ;
    dPres = ( LOG10(MaxPres) - LOG10(MinPres) ) / nPres ;

    fseek(fp, 8, SEEK_CUR); // skip blank
    Table_LogRho = ( double * ) malloc( sizeOfTable );
    if ( Table_LogRho == NULL ) Aux_Error( ERROR_INFO, "Could not allocate memory for Table_LogRho\n");
    fread( Table_LogRho, sizeOfTable, 1, fp );

    fseek(fp, 8, SEEK_CUR); // skip blank
    Table_LogEner = ( double *) malloc( sizeOfTable );
    if ( Table_LogEner == NULL ) Aux_Error( ERROR_INFO, "Could not allocate memory for Table_LogEner\n");
    fread( Table_LogEner, sizeOfTable, 1, fp );

    fseek(fp, 8, SEEK_CUR); // skip blank
    Table_DensEint2Temp = ( double * ) malloc( sizeOfTable );
    if ( Table_DensEint2Temp == NULL ) Aux_Error( ERROR_INFO, "Could not allocate memory for Table_DensEint2Temp\n");
    fread( Table_DensEint2Temp, sizeOfTable, 1, fp );

    fseek(fp, 8, SEEK_CUR); // skip blank
    Table_DensEint2Pres = ( double * ) malloc( sizeOfTable );
    if ( Table_DensEint2Pres == NULL ) Aux_Error( ERROR_INFO, "Could not allocate memory for Table_DensEint2Pres\n");
    fread( Table_DensEint2Pres, sizeOfTable, 1, fp );

    fseek(fp, 8, SEEK_CUR); // skip blank
    fseek(fp, sizeOfTable, SEEK_CUR); // skip S_EoS_Table, because S_EoS_Table is useless.

    fseek(fp, 8, SEEK_CUR); // skip blank
    Table_DensEint2Cs = ( double * ) malloc( sizeOfTable );
    if ( Table_DensEint2Cs == NULL ) Aux_Error( ERROR_INFO, "Could not allocate memory for Table_DensEint2Cs\n");
    fread( Table_DensEint2Cs, sizeOfTable, 1, fp );

    fseek(fp, 8, SEEK_CUR); // skip blank
    fseek(fp, sizeOfTable, SEEK_CUR); // skip xH_EoS, because xH_EoS is useless.

    fseek(fp, 8, SEEK_CUR); // skip blank
    fseek(fp, sizeOfTable, SEEK_CUR); // skip xH2_EoS, because xH2_EoS is useless.

    fseek(fp, 8, SEEK_CUR); // skip blank
    fseek(fp, sizeOfTable, SEEK_CUR); // skip xHe_EoS, because xHe_EoS is useless.

    fseek(fp, 8, SEEK_CUR); // skip blank
    fseek(fp, sizeOfTable, SEEK_CUR); // skip xHep_EoS, because xHep_EoS is useless.

    fclose(fp);

    for ( int i = 0; i < nRho * nEner; i++ )
    {
        Table_LogRho[i] = LOG10( Table_LogRho[i] );
        Table_LogEner[i] = LOG10( Table_LogEner[i] );
    }
    
    h_EoS_Table[SAUMON_CHABRIER_TAB_LOG_RHO] = Table_LogRho;
    h_EoS_Table[SAUMON_CHABRIER_TAB_LOG_ENER] = Table_LogEner;
    h_EoS_Table[SAUMON_CHABRIER_TAB_DENSEINT2TEMP] = Table_DensEint2Temp;
    h_EoS_Table[SAUMON_CHABRIER_TAB_DENSEINT2PRES] = Table_DensEint2Pres;
    h_EoS_Table[SAUMON_CHABRIER_TAB_DENSEINT2CS] = Table_DensEint2Cs;

    if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

}
#endif // if ( MODEL == HYDRO && EOS == EOS_BAROTROP1D )