#include <cstdio>

#define SQR(  a )       ( (a)*(a)     )
#define CUBE( a )       ( (a)*(a)*(a) )

int main()
{
    const double n_0        =  100;      // intial density (cc)
    const double v_flow     =    2;      // the magnitude of flow velocity (km/s)
    const int NX            =   48;      // number of base-level cells along x
    const int NY            =   48;      // number of base-level cells along y
    const int NZ            =  592;      // number of base-level cells along z
    const double GAMMA      =  5/3;
    const int NFIELD        =    5;
    
    const double UNIT_L     =    3.08567758149e18;    // pc to cm
    const double UNIT_V     =                1000;    // km/s to cm/s
    const double UNIT_D     =      3.84703042e-24;    // 1/cm3 to g/cm3
    const double UNIT_M     = UNIT_D*CUBE(UNIT_L);
    const double UNIT_E     =  UNIT_M*SQR(UNIT_V);
    
    const double MOLECULAR_WEIGHT =           2.3;
    const double ISO_TEMP         =            10;    // Isothermal temperature (K)
    const double Const_kB        = 1.38064852e-16;    // Boltzmann constant in erg/K
    const double Const_amu      = 1.660539040e-24;    // atomic mass unit, proton mass in g

    float (*IC)[NZ][NY][NX] = new float [NFIELD][NZ][NY][NX];

    const double rho_0      = Const_amu*MOLECULAR_WEIGHT*n_0;   // convert to mass density (g/cm3)
    const double Cs2        = ( Const_kB*ISO_TEMP/UNIT_E ) / ( MOLECULAR_WEIGHT*Const_amu/UNIT_M );
    const double Pres       = Cs2*rho_0;
    const double Eint       = (double)1.0e4*Pres;
    
    double MomX, MomY, MomZ, Etot;
    for (int k=0; k<NZ; k++)
    for (int j=0; j<NY; j++)
    for (int i=0; i<NX; i++)
    {
        if (j <= NY/2)
           MomX = rho_0*v_flow;
        else
           MomX = -rho_0*v_flow;
        MomY = 0;
        MomZ = 0;

        Etot = (double)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) ) / rho_0;
        Etot += Eint;

        IC[0][k][j][i] = n_0;    // mass density
        IC[1][k][j][i] = MomX;   // momentum density x
        IC[2][k][j][i] = MomY;   // momentum density y
        IC[3][k][j][i] = MomZ;   // momentum density z
        IC[4][k][j][i] = Etot;   // total energy density
    }

    FILE *File = fopen( "UM_IC", "wb" );
    fwrite( IC, sizeof(float), NFIELD*NZ*NY*NX, File );
    fclose( File );

    delete [] IC;
}

