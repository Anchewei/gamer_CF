#include "EoS_Saumon_Chabrier.h"
#include "GAMER.h"

#if (MODEL == HYDRO && EOS == EOS_SAUMON_CHABRIER)

extern int nRho;
extern int nEner;

void Fill_Hole(double *Table, const int Table_Index)
{
    if (MPI_Rank == 0)
        Aux_Message(stdout, "%s ... \n", __FUNCTION__);

    if (nRho != 625)
        Aux_Error(ERROR_INFO, "Number of Rho = %d\n", nRho);
    if (nEner != 1000)
        Aux_Error(ERROR_INFO, "Number of Ener = %d\n", nEner);

    switch (Table_Index)
    {
    case SAUMON_CHABRIER_TAB_DENSEINT2TEMP:
        if (MPI_Rank == 0)
            Aux_Message(stdout, "   Fill Table Temp ... \n");
        break;

    case SAUMON_CHABRIER_TAB_DENSEINT2PRES:
        if (MPI_Rank == 0)
            Aux_Message(stdout, "   Fill Table Pres ... \n");
        break;
    case SAUMON_CHABRIER_TAB_DENSEINT2CS:
        if (MPI_Rank == 0)
            Aux_Message(stdout, "   Fill Table Cs ... \n");
        break;

    default:
        Aux_Error(ERROR_INFO, "Not supported fill table index = %d", Table_Index);
        break;
    }

    int ii = 0, jj = 0, kk = 0, hh = 0, ee = 0, gg = 0;

    for (int k = 0; k < 5; k++)
    {
        ii = 0;
        jj = 0;
        kk = 0;
        hh = 0;
        ee = 0;
        gg = 0;
        double ww, xx, yy, zz;

        for (int ir = 1; ir < nRho - 1; ir++)
        {
            for (int ie = 1; ie < nEner - 1; ie++)
            {
                if (Table[ie * nRho + ir] == 0.0)
                {
                    ii++;
                    xx = Table[(ie + 1) * nRho + ir] * Table[(ie - 1) * nRho + ir] *
                         Table[ie * nRho + ir - 1] * Table[ie * nRho + ir + 1];
                    yy = Table[(ie + 1) * nRho + ir + 1] * Table[(ie - 1) * nRho + ir + 1] *
                         Table[(ie - 1) * nRho + ir - 1] * Table[(ie + 1) * nRho + ir - 1];

                    if (ie > 1 && ie < nEner - 2 && ir > 1 && ir < nRho - 2)
                    {
                        ww = Table[(ie + 2) * nRho + ir] * Table[(ie - 2) * nRho + ir] *
                             Table[ie * nRho + ir - 2] * Table[ie * nRho + ir + 2];
                    }
                    else
                        ww = 0.0;

                    if (ie > 2 && ie < nEner - 3 && ir > 2 && ir < nRho - 3)
                    {
                        zz = Table[(ie + 3) * nRho + ir + 3] * Table[(ie - 3) * nRho + ir - 3] *
                             Table[(ie + 3) * nRho + ir - 3] * Table[(ie - 3) * nRho + ir + 3];
                    }
                    else
                        zz = 0.0;

                    if (xx != 0.0)
                    {
                        jj++;
                        Table[ie * nRho + ir] = 0.25 * (Table[(ie + 1) * nRho + ir] + Table[(ie - 1) * nRho + ir] +
                                                        Table[ie * nRho + ir - 1] + Table[ie * nRho + ir + 1]);
                    }
                    else if (yy != 0.0 && k >= 0)
                    {
                        kk++;
                        Table[ie * nRho + ir] = 0.25 * (Table[(ie + 1) * nRho + ir + 1] + Table[(ie - 1) * nRho + ir + 1] +
                                                        Table[(ie + 1) * nRho + ir - 1] + Table[(ie - 1) * nRho + ir - 1]);
                    }
                    else if (ww != 0 && k > 0)
                    {
                        ee++;
                        Table[ie * nRho + ir] = 0.25 * (Table[(ie + 2) * nRho + ir] + Table[(ie - 2) * nRho + ir] +
                                                        Table[ie * nRho + ir - 2] + Table[ie * nRho + ir + 2]);
                    }
                    else if (zz != 0 && k > 1)
                    {
                        hh++;
                        Table[ie * nRho + ir] = 0.25 * (Table[(ie + 3) * nRho + ir + 3] + Table[(ie - 3) * nRho + ir + 3] +
                                                        Table[(ie + 3) * nRho + ir - 3] + Table[(ie - 3) * nRho + ir - 3]);
                    }
                    else
                        gg++;
                }
            }
        }
        if (MPI_Rank == 0)
            Aux_Message(stdout, "   ii = %5d, jj = %5d, kk = %5d, ee = %5d, hh = %5d, gg = %5d, iter = %5d\n", ii, jj, kk, ee, hh, gg, k);
    }

    switch (Table_Index)
    {
    case SAUMON_CHABRIER_TAB_DENSEINT2TEMP:
        if (MPI_Rank == 0)
            Aux_Message(stdout, "   Fill Table Temp ... done\n");
        break;

    case SAUMON_CHABRIER_TAB_DENSEINT2PRES:
        if (MPI_Rank == 0)
            Aux_Message(stdout, "   Fill Table Pres ... done\n");
        break;
    case SAUMON_CHABRIER_TAB_DENSEINT2CS:
        if (MPI_Rank == 0)
            Aux_Message(stdout, "   Fill Table Cs ... done\n");
        break;

    default:
        Aux_Error(ERROR_INFO, "Not supported fill table index = %d", Table_Index);
        break;
    }

    if (MPI_Rank == 0)
        Aux_Message(stdout, "%s ... done\n", __FUNCTION__);
}
#endif