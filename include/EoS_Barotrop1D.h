#ifndef __EOS_BAROTROP1D_H__
#define __EOS_BAROTROP1D_H__

#include "CUFLU.h"
#define BAROTROP1D_TABLE_NVAR
#define BAROTROP1D_TABLE_NPTR

// auxiliary array indices
#define BAROTROP1D_AUX_NRHO         0 // AuxArray_Int: Number of row in Rho Table

#define BAROTROP1D_AUX_MAXRHO       0 // AuxArray_Flt: Max Rho
#define BAROTROP1D_AUX_MINRHO       1 // AuxArray_Flt: Min Rho
#define BAROTROP1D_AUX_DRHO         2 // AuxArray_Flt: dRho
#define BAROTROP1D_AUX_GAMMA        3 // AuxArray_Flt: GAMMA
#define BAROTROP1D_AUX_GAMMA_M1     4 // AuxArray_Flt: GAMMA - 1.0
#define BAROTROP1D_AUX__GAMMA_M1    5 // AuxArray_Flt: 1 / ( GAMMA - 1.0 );
#define BAROTROP1D_AUX_DENS2CGS     6 // AuxArray_Flt: convert density to cgs
#define BAROTROP1D_AUX_M_KB         7 // AuxArray_Flt: mean molecular weight*atomic mass unit/ Bolzmann constant*(UNIT_E/UNIT_M)
#define BAROTROP1D_AUX__M_KB        8 // AuxArray_Flt: 1 / ( mean molecular weight*atomic mass unit/ Bolzmann constant*(UNIT_E/UNIT_M) )

// table indices
#define BAROTROP1D_TAB_ALL          0 // All Tables
#define BAROTROP1D_TAB_RHO          1 // Rho Tables
#define BAROTROP1D_TAB_TEMP         2 // Temp Tables

// variable indices in the Barotrop1D EoS table and auxiliary table

#endif