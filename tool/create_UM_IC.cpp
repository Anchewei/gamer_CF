#include <cstdio>

int main()
{
    const double n_0        =  100;      // intial density (cc)
    const double v_flow     =    2;      // the magnitude of flow velocity (km/s)
    const int NX            =   64;      //number of base-level cells along x
    const int NY            =   64;      //number of base-level cells along y
    const int NZ            =  608;      //number of base-level cells along z
    const int NFIELD        =    5;

    float (*IC)[NZ][NY][NX] = new float [NFIELD][NZ][NY][NX];

    const double rho_0      = 3.84703042e-24*n_0;   // convert to mass density (g/cm3), assuming mu_H = 2.3
    
    for (int k=0; k<NZ; k++)
    for (int j=0; j<NY; j++)
    for (int i=0; i<NX; i++)
    {
       IC[0][k][j][i] = n_0;   // mass density
       if (j <= NY/2)
           IC[1][k][j][i] = rho_0*v_flow;   // momentum density x
       else
           IC[1][k][j][i] = -rho_0*v_flow;   // momentum density x
       IC[2][k][j][i] = 0.0;   // momentum density y
       IC[3][k][j][i] = 0.0;   // momentum density z
       IC[4][k][j][i] = 2.0;   // total energy density
    }

    FILE *File = fopen( "UM_IC", "wb" );
    fwrite( IC, sizeof(float), NFIELD*NZ*NY*NX, File );
    fclose( File );

    delete [] IC;
}

