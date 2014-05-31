// written by Szymon Winczewski
// based on LAMMPS implementation of AIREBO force field

#include "airebo_force_field.h"

using namespace std;


double AIREBOForceField::gSpline(double costh, double Nij, int typei)
{
    double coeffs[6], g1, g2, cut, g;
    int i, j;

    i = 0;
    j = 0;
    g = 0.0;
    cut = 0.0;

// !!!komentarz: central atom is Carbon
    if ( typei == 0 )
    {
        if ( costh < gCdom[0] )
            costh = gCdom[0];
        if ( costh > gCdom[4] )
            costh = gCdom[4];
        if ( Nij >= NCmax )
        {
            for (i = 0; i < 4; i++)
            {
                if ( ( costh >= gCdom[i] ) && ( costh <= gCdom[i+1] ) )
                {
                    for (j = 0; j < 6; j++)
                        coeffs[j] = gC2[i][j];
                }
            }
            g2 = Sp5th(costh, coeffs);
            g = g2;
        }
        if ( Nij <= NCmin )
        {
            for (i = 0; i < 4; i++)
            {
                if ( ( costh >= gCdom[i] ) && ( costh <= gCdom[i+1] ) )
                {
                    for (j = 0; j < 6; j++)
                        coeffs[j] = gC1[i][j];
                }
            }
            g1 = Sp5th(costh, coeffs);
            g = g1;
        }
        if ( ( Nij > NCmin ) && ( Nij < NCmax ) )
        {
            for (i = 0; i < 4; i++)
            {
                if ( ( costh >= gCdom[i] ) && ( costh <= gCdom[i+1] ) )
                {
                    for (j = 0; j < 6; j++)
                        coeffs[j] = gC1[i][j];
                }
            }
            g1 = Sp5th(costh, coeffs);
            for (i = 0; i < 4; i++)
            {
                if ( ( costh >= gCdom[i] ) && ( costh <= gCdom[i+1] ) )
                {
                    for (j = 0; j < 6; j++)
                        coeffs[j] = gC2[i][j];
                }
            }
            g2 = Sp5th(costh, coeffs);
            cut = Sp(Nij, NCmin, NCmax);
            g = g2 + cut * ( g1 - g2 );
        }
    }

// !!!komentarz: central atom is Hydrogen
    if ( typei == 1 )
    {
        if ( costh < gHdom[0] )
            costh = gHdom[0];
        if ( costh > gHdom[3] )
            costh = gHdom[3];
        for (i = 0; i < 3; i++)
        {
            if ( (costh >= gHdom[i] ) && ( costh <= gHdom[i+1] ) )
            {
                for (j = 0; j < 6; j++)
                    coeffs[j] = gH[i][j];
            }
        }
        g = Sp5th(costh, coeffs);
    }

    return g;
}


double AIREBOForceField::PijSpline(double NijC, double NijH, int typei, int typej)
{
    int x, y, i;
    double Pij, coeffs[16];

    x = 0;
    y = 0;

// !!!komentarz: if the inputs are out of bounds set them back to a point in bounds
    if ( ( typei == 0 ) && ( typej == 0 ) )
    {
        if ( NijC < pCCdom[0][0] )
            NijC = pCCdom[0][0];
        if ( NijC > pCCdom[0][1] )
            NijC = pCCdom[0][1];
        if ( NijH < pCCdom[1][0] )
            NijH = pCCdom[1][0];
        if ( NijH > pCCdom[1][1] )
            NijH = pCCdom[1][1];

        if ( ( fabs( NijC - floor(NijC) ) < TOL ) && ( fabs( NijH - floor(NijH) ) < TOL ) )
        {
            Pij = PCCf[(int) NijC][(int) NijH];
            return Pij;
        }

        x = (int) (floor(NijC));
        y = (int) (floor(NijH));
        for (i = 0; i < 16; i++)
            coeffs[i] = pCC[x][y][i];
        Pij = Spbicubic(NijC, NijH, coeffs);
        return Pij;
    }

// !!!komentarz: if the inputs are out of bounds set them back to a point in bounds
    if ( ( typei == 0 ) && ( typej == 1 ) )
    {
        if ( NijC < pCHdom[0][0] )
            NijC=pCHdom[0][0];
        if ( NijC > pCHdom[0][1] )
            NijC=pCHdom[0][1];
        if ( NijH < pCHdom[1][0] )
            NijH=pCHdom[1][0];
        if ( NijH > pCHdom[1][1] )
            NijH=pCHdom[1][1];

        if ( ( fabs( NijC - floor(NijC) ) < TOL ) && ( fabs( NijH - floor(NijH) ) < TOL ) )
        {
            Pij = PCHf[(int) NijC][(int) NijH];
            return Pij;
        }

        x = (int) (floor(NijC));
        y = (int) (floor(NijH));
        for (i = 0; i < 16; i++)
            coeffs[i] = pCH[x][y][i];
        Pij = Spbicubic(NijC, NijH, coeffs);
        return Pij;
    }


    if ( ( typei == 1 ) && ( typej == 0 ) )
        return 0.0;
    if ( ( typei == 1 ) && ( typej == 1 ) )
        return 0.0;

    return 0.0;
}


double AIREBOForceField::piRCSpline(double Nij, double Nji, double Nijconj, int typei, int typej)
{
    int x, y, z, i, done;
    double piRC, coeffs[64];

    x = 0;
    y = 0;
    z = 0;
    i = 0;
    done = 0;
    piRC = 0.0;

    for (i = 0; i < 64; i++)
        coeffs[i] = 0.0;

    if ( ( typei == 0 ) && ( typej == 0 ) )
    {
// !!!komentarz: if the inputs are out of bounds set them back to a point in bounds
        if ( Nij < piCCdom[0][0] )
            Nij = piCCdom[0][0];
        if ( Nij > piCCdom[0][1] )
            Nij = piCCdom[0][1];
        if ( Nji < piCCdom[1][0] )
            Nji = piCCdom[1][0];
        if ( Nji > piCCdom[1][1] )
            Nji = piCCdom[1][1];
        if ( Nijconj < piCCdom[2][0] )
            Nijconj = piCCdom[2][0];
        if ( Nijconj > piCCdom[2][1] )
            Nijconj = piCCdom[2][1];

        if ( ( fabs( Nij - floor(Nij) ) < TOL ) && ( fabs( Nji - floor(Nji) ) < TOL ) &&
            ( fabs( Nijconj - floor(Nijconj) ) < TOL ) )
        {
            piRC = piCCf[(int) Nij][(int) Nji][(int) Nijconj];
            done = 1;
        }

        if ( done == 0 )
        {
// !!!uwaga: ponizej niejednoznaczne warunki
// if ( Nij >= (double) i && Nij <= (double) i+1 || Nij == (double) i )
// if ( Nji >= (double) i && Nji <= (double) i+1 || Nji == (double) i )
// if ( Nijconj >= (double) i && Nijconj <= (double) i+1 || Nijconj == (double) i )
            for (i = 0; i < piCCdom[0][1]; i++)
                if ( ( ( Nij >= (double) i ) && ( Nij <= (double) i+1 ) ) || ( Nij == (double) i ) )
                    x = i;
            for (i = 0; i < piCCdom[1][1]; i++)
                if ( ( ( Nji >= (double) i ) && ( Nji <= (double) i+1 ) ) || ( Nji == (double) i ) )
                    y = i;
            for (i = 0; i < piCCdom[2][1]; i++)
                if ( ( ( Nijconj >= (double) i ) && ( Nijconj <= (double) i+1 ) ) || ( Nijconj == (double) i ) )
                    z = i;
            for (i = 0; i < 64; i++)
                coeffs[i] = piCC[x][y][z][i];
            piRC = Sptricubic(Nij, Nji, Nijconj, coeffs);
        }
    }

// !!!komentarz: CH interaction
    if ( ( ( typei == 0 ) && ( typej == 1 ) ) || ( ( typei == 1 ) && ( typej == 0 ) ) )
    {
// !!!komentarz: if the inputs are out of bounds set them back to a point in bounds
        if ( ( Nij < piCHdom[0][0] ) || ( Nij > piCHdom[0][1] ) ||
            ( Nji < piCHdom[1][0] ) || ( Nji > piCHdom[1][1] ) ||
            ( Nijconj < piCHdom[2][0] ) || ( Nijconj > piCHdom[2][1] ) )
        {
            if ( Nij < piCHdom[0][0] )
                Nij = piCHdom[0][0];
            if ( Nij > piCHdom[0][1] )
                Nij = piCHdom[0][1];
            if ( Nji < piCHdom[1][0] )
                Nji = piCHdom[1][0];
            if ( Nji > piCHdom[1][1] )
                Nji = piCHdom[1][1];
            if ( Nijconj < piCHdom[2][0] )
                Nijconj = piCHdom[2][0];
            if ( Nijconj > piCHdom[2][1] )
                Nijconj = piCHdom[2][1];
        }

        if ( ( fabs( Nij - floor(Nij) ) < TOL ) && ( fabs( Nji - floor(Nji) ) < TOL ) &&
            ( fabs( Nijconj - floor(Nijconj) ) < TOL ) )
        {
            piRC = piCHf[(int) Nij][(int) Nji][(int) Nijconj];
            done = 1;
        }

        if ( done == 0 )
        {
            for (i = 0; i < piCHdom[0][1]; i++)
                if ( Nij >= i && Nij <= i+1 )
                    x = i;
            for (i = 0; i < piCHdom[1][1]; i++)
                if ( Nji >= i && Nji <= i+1 )
                    y = i;
            for (i = 0; i < piCHdom[2][1]; i++)
                if ( Nijconj >= i && Nijconj <= i+1)
                    z = i;
            for (i = 0; i < 64; i++)
                coeffs[i] = piCH[x][y][z][i];
            piRC = Sptricubic(Nij, Nji, Nijconj, coeffs);
        }
    }

    if ( ( typei == 1 ) && ( typej == 1 ) )
    {
        if ( ( Nij < piHHdom[0][0] ) || ( Nij > piHHdom[0][1] ) ||
            ( Nji < piHHdom[1][0] ) || ( Nji > piHHdom[1][1] ) ||
            ( Nijconj < piHHdom[2][0] ) || ( Nijconj > piHHdom[2][1] ) )
        {
            Nij = 0.0;
            Nji = 0.0;
            Nijconj = 0.0;
        }
        if ( ( fabs( Nij - floor(Nij) ) < TOL ) && ( fabs( Nji - floor(Nji) ) < TOL ) &&
            ( fabs( Nijconj - floor(Nijconj) ) < TOL ) )
        {
            piRC = piHHf[(int) Nij][(int) Nji][(int) Nijconj];
            done = 1;
        }
        if ( done == 0 )
        {
            for (i = 0; i < piHHdom[0][1]; i++)
                if ( Nij >=i && Nij <= i+1 )
                    x = i;
            for (i = 0; i < piHHdom[1][1]; i++)
                if ( Nji >=i && Nji <= i+1 )
                    y = i;
            for (i = 0; i < piHHdom[2][1]; i++)
                if ( Nijconj >=i && Nijconj <= i+1 )
                    z = i;
            for (i = 0; i < 64; i++)
                coeffs[i] = piHH[x][y][z][i];
            piRC = Sptricubic(Nij, Nji, Nijconj, coeffs);
        }
    }

    return piRC;
}


double AIREBOForceField::TijSpline(double Nij, double Nji, double Nijconj)
{
    int x, y, z, i, done;
    double Tijf, coeffs[64];

    x = 0;
    y = 0;
    z = 0;
    i = 0;
    Tijf = 0.0;
    done = 0;
    for (i = 0; i < 64; i++)
        coeffs[i] = 0.0;

// !!!komentarz: if the inputs are out of bounds set them back to a point in bounds
    if ( Nij < Tijdom[0][0] )
        Nij=Tijdom[0][0];
    if ( Nij > Tijdom[0][1] )
        Nij=Tijdom[0][1];
    if ( Nji < Tijdom[1][0] )
        Nji=Tijdom[1][0];
    if ( Nji > Tijdom[1][1] )
        Nji=Tijdom[1][1];
    if ( Nijconj < Tijdom[2][0] )
        Nijconj=Tijdom[2][0];
    if ( Nijconj > Tijdom[2][1] )
        Nijconj=Tijdom[2][1];

    if ( ( fabs( Nij-floor(Nij) ) < TOL ) && ( fabs( Nji - floor(Nji) ) < TOL ) &&
        ( fabs( Nijconj - floor(Nijconj) ) < TOL ) )
    {
        Tijf = Tf[(int) Nij][(int) Nji][(int) Nijconj];
        done = 1;
    }

    if ( done == 0 )
    {
        for (i = 0; i < Tijdom[0][1]; i++)
            if ( Nij >= i && Nij <= i+1)
                x = i;
        for (i = 0; i < Tijdom[1][1]; i++)
            if ( Nji >= i && Nji <= i+1)
                y = i;
        for (i = 0; i < Tijdom[2][1]; i++)
            if ( Nijconj >= i && Nijconj <= i+1)
                z = i;
        for (i = 0; i < 64; i++)
            coeffs[i] = Tijc[x][y][z][i];
        Tijf = Sptricubic(Nij, Nji, Nijconj, coeffs);
    }

    return Tijf;
}


double AIREBOForceField::Sp5th(double x, double *coeffs)
{
    double f;
    const double x2 = x * x;
    const double x3 = x2 * x;

    f  = coeffs[0];
    f += coeffs[1] * x;
    f += coeffs[2] * x2;
    f += coeffs[3] * x3;
    f += coeffs[4] * x2 * x2;
    f += coeffs[5] * x2 * x3;

    return f;
}


double AIREBOForceField::Spbicubic(double x, double y, double *coeffs)
{
    double f, xn, yn, c;
    int i, j;

    f = 0.0;

    xn = 1.0;
    for (i = 0; i < 4; i++)
    {
        yn = 1.0;
        for (j = 0; j < 4; j++)
        {
            c = coeffs[i * 4 + j];
            f += c * xn * yn;
            yn *= y;
        }
        xn *= x;
    }

    return f;
}


double AIREBOForceField::Sptricubic(double x, double y, double z, double *coeffs)
{
    double f, xn, yn, zn, c;
    int i, j, k;

    f = 0.0;

    xn = 1.0;
    for (i = 0; i < 4; i++)
    {
        yn = 1.0;
        for (j = 0; j < 4; j++)
        {
            zn = 1.0;
            for (k = 0; k < 4; k++)
            {
                c = coeffs[16 * i + 4 * j + k];
                f += c * xn * yn * zn;
                zn *= z;
            }
            yn *= y;
        }
        xn *= x;
    }

    return f;
}
