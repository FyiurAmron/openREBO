/*
    File:   OpenREBO_splines.cpp
    Author: vaxquis
    based on LAMMPS implementation of AIREBO potential
 */

#include "OpenREBO.h"

namespace OpenREBO {
  const int DIM = 4;

  double AIREBO::gSpline( double costh, double Nij, int typei ) {
    double coeffs[6], g1, g2;
    double cut = 0.0, g = 0.0;

    int i, j;

    if ( typei == 0 ) { // central atom is Carbon
      if ( costh < gCdom[0] )
        costh = gCdom[0];
      if ( costh > gCdom[4] )
        costh = gCdom[4];
      if ( Nij >= NCmax ) {
        for( i = 0; i < 4; i++ ) {
          if ( ( costh >= gCdom[i] ) && ( costh <= gCdom[i + 1] ) ) {
            for( j = 0; j < 6; j++ )
              coeffs[j] = gC2[i][j];
          }
        }
        g2 = Sp5th( costh, coeffs );
        g = g2;
      }
      if ( Nij <= NCmin ) {
        for( i = 0; i < 4; i++ ) {
          if ( ( costh >= gCdom[i] ) && ( costh <= gCdom[i + 1] ) ) {
            for( j = 0; j < 6; j++ )
              coeffs[j] = gC1[i][j];
          }
        }
        g1 = Sp5th( costh, coeffs );
        g = g1;
      }
      if ( ( Nij > NCmin ) && ( Nij < NCmax ) ) {
        for( i = 0; i < 4; i++ ) {
          if ( ( costh >= gCdom[i] ) && ( costh <= gCdom[i + 1] ) ) {
            for( j = 0; j < 6; j++ )
              coeffs[j] = gC1[i][j];
          }
        }
        g1 = Sp5th( costh, coeffs );
        for( i = 0; i < 4; i++ ) {
          if ( ( costh >= gCdom[i] ) && ( costh <= gCdom[i + 1] ) ) {
            for( j = 0; j < 6; j++ )
              coeffs[j] = gC2[i][j];
          }
        }
        g2 = Sp5th( costh, coeffs );
        cut = SpNC( Nij );
        g = g2 + cut * ( g1 - g2 );
      }
    } else /*if ( typei == 1 )*/ { // central atom is Hydrogen
      if ( costh < gHdom[0] )
        costh = gHdom[0];
      if ( costh > gHdom[3] )
        costh = gHdom[3];
      for( i = 0; i < 3; i++ ) {
        if ( ( costh >= gHdom[i] ) && ( costh <= gHdom[i + 1] ) ) {
          for( j = 0; j < 6; j++ )
            coeffs[j] = gH[i][j];
        }
      }
      g = Sp5th( costh, coeffs );
    }

    return g;
  }

  double AIREBO::PijSpline( double NijC, double NijH, int typei, int typej ) {
    int x, y, i;
    double Pij, coeffs[16];

    x = 0;
    y = 0;

    if ( typei == 1 )
      return 0.0;

    // if the inputs are out of bounds set them back to a point in bounds
    if ( typej == 0 ) {
      if ( NijC < pCCdom[0][0] )
        NijC = pCCdom[0][0];
      if ( NijC > pCCdom[0][1] )
        NijC = pCCdom[0][1];
      if ( NijH < pCCdom[1][0] )
        NijH = pCCdom[1][0];
      if ( NijH > pCCdom[1][1] )
        NijH = pCCdom[1][1];

      if ( ( fabs( NijC - floor( NijC ) ) < TOL ) && ( fabs( NijH - floor( NijH ) ) < TOL ) ) {
        Pij = PCCf[(int) NijC][(int) NijH];
        return Pij;
      }

      x = (int) ( floor( NijC ) );
      y = (int) ( floor( NijH ) );
      for( i = 0; i < 16; i++ )
        coeffs[i] = pCC[x][y][i];
      Pij = Spbicubic( NijC, NijH, coeffs );
      return Pij;
    }

    // else ( typej == 1 )
    if ( NijC < pCHdom[0][0] )
      NijC = pCHdom[0][0];
    if ( NijC > pCHdom[0][1] )
      NijC = pCHdom[0][1];
    if ( NijH < pCHdom[1][0] )
      NijH = pCHdom[1][0];
    if ( NijH > pCHdom[1][1] )
      NijH = pCHdom[1][1];

    if ( ( fabs( NijC - floor( NijC ) ) < TOL ) && ( fabs( NijH - floor( NijH ) ) < TOL ) ) {
      Pij = PCHf[(int) NijC][(int) NijH];
      return Pij;
    }

    x = (int) ( floor( NijC ) );
    y = (int) ( floor( NijH ) );
    for( i = 0; i < 16; i++ )
      coeffs[i] = pCH[x][y][i];
    Pij = Spbicubic( NijC, NijH, coeffs );
    return Pij;
  }

  double AIREBO::piRCSpline( double Nij, double Nji, double Nijconj, int typei, int typej ) {
    int x, y, z, i;
    double piRC, coeffs[64] = { };

    x = 0;
    y = 0;
    z = 0;
    piRC = 0.0;

    if ( ( typei == 0 ) && ( typej == 0 ) ) { // CC
      // if the inputs are out of bounds set them back to a point in bounds
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

      if ( ( fabs( Nij - floor( Nij ) ) < TOL ) && ( fabs( Nji - floor( Nji ) ) < TOL ) &&
              ( fabs( Nijconj - floor( Nijconj ) ) < TOL ) )
        return piCCf[(int) Nij][(int) Nji][(int) Nijconj];

      // !!!uwaga: ponizej niejednoznaczne warunki (?)
      // if ( Nij >= (double) i && Nij <= (double) i+1 || Nij == (double) i )
      // if ( Nji >= (double) i && Nji <= (double) i+1 || Nji == (double) i )
      // if ( Nijconj >= (double) i && Nijconj <= (double) i+1 || Nijconj == (double) i )
      for( i = 0; i < piCCdom[0][1]; i++ )
        if ( ( ( Nij >= (double) i ) && ( Nij <= (double) i + 1 ) ) || ( Nij == (double) i ) )
          x = i;
      for( i = 0; i < piCCdom[1][1]; i++ )
        if ( ( ( Nji >= (double) i ) && ( Nji <= (double) i + 1 ) ) || ( Nji == (double) i ) )
          y = i;
      for( i = 0; i < piCCdom[2][1]; i++ )
        if ( ( ( Nijconj >= (double) i ) && ( Nijconj <= (double) i + 1 ) ) || ( Nijconj == (double) i ) )
          z = i;
      for( i = 0; i < 64; i++ )
        coeffs[i] = piCC[x][y][z][i];
      return Sptricubic( Nij, Nji, Nijconj, coeffs );
    } else if ( ( typei == 1 ) && ( typej == 1 ) ) { // HH
      if ( ( Nij < piHHdom[0][0] ) || ( Nij > piHHdom[0][1] ) ||
              ( Nji < piHHdom[1][0] ) || ( Nji > piHHdom[1][1] ) ||
              ( Nijconj < piHHdom[2][0] ) || ( Nijconj > piHHdom[2][1] ) ) {
        Nij = 0.0;
        Nji = 0.0;
        Nijconj = 0.0;
      }
      if ( ( fabs( Nij - floor( Nij ) ) < TOL ) && ( fabs( Nji - floor( Nji ) ) < TOL ) &&
              ( fabs( Nijconj - floor( Nijconj ) ) < TOL ) )
        return piHHf[(int) Nij][(int) Nji][(int) Nijconj];
      for( i = 0; i < piHHdom[0][1]; i++ )
        if ( Nij >= i && Nij <= i + 1 )
          x = i;
      for( i = 0; i < piHHdom[1][1]; i++ )
        if ( Nji >= i && Nji <= i + 1 )
          y = i;
      for( i = 0; i < piHHdom[2][1]; i++ )
        if ( Nijconj >= i && Nijconj <= i + 1 )
          z = i;
      for( i = 0; i < 64; i++ )
        coeffs[i] = piHH[x][y][z][i];
      return Sptricubic( Nij, Nji, Nijconj, coeffs );
    }
    // else CH
    // if the inputs are out of bounds set them back to a point in bounds
    if ( ( Nij < piCHdom[0][0] ) || ( Nij > piCHdom[0][1] ) ||
            ( Nji < piCHdom[1][0] ) || ( Nji > piCHdom[1][1] ) ||
            ( Nijconj < piCHdom[2][0] ) || ( Nijconj > piCHdom[2][1] ) ) {
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

    if ( ( fabs( Nij - floor( Nij ) ) < TOL ) && ( fabs( Nji - floor( Nji ) ) < TOL ) &&
            ( fabs( Nijconj - floor( Nijconj ) ) < TOL ) )
      return piCHf[(int) Nij][(int) Nji][(int) Nijconj];

    for( i = 0; i < piCHdom[0][1]; i++ )
      if ( Nij >= i && Nij <= i + 1 )
        x = i;
    for( i = 0; i < piCHdom[1][1]; i++ )
      if ( Nji >= i && Nji <= i + 1 )
        y = i;
    for( i = 0; i < piCHdom[2][1]; i++ )
      if ( Nijconj >= i && Nijconj <= i + 1 )
        z = i;
    for( i = 0; i < 64; i++ )
      coeffs[i] = piCH[x][y][z][i];
    return Sptricubic( Nij, Nji, Nijconj, coeffs );
  }

  double AIREBO::TijSpline( double Nij, double Nji, double Nijconj ) {
    bool done;
    int x, y, z, i;
    double Tijf, coeffs[64] = { };

    x = 0;
    y = 0;
    z = 0;
    i = 0;
    Tijf = 0.0;
    done = 0;

    // if the inputs are out of bounds set them back to a point in bounds
    if ( Nij < Tijdom[0][0] )
      Nij = Tijdom[0][0];
    if ( Nij > Tijdom[0][1] )
      Nij = Tijdom[0][1];
    if ( Nji < Tijdom[1][0] )
      Nji = Tijdom[1][0];
    if ( Nji > Tijdom[1][1] )
      Nji = Tijdom[1][1];
    if ( Nijconj < Tijdom[2][0] )
      Nijconj = Tijdom[2][0];
    if ( Nijconj > Tijdom[2][1] )
      Nijconj = Tijdom[2][1];

    if ( ( fabs( Nij - floor( Nij ) ) < TOL ) && ( fabs( Nji - floor( Nji ) ) < TOL ) &&
            ( fabs( Nijconj - floor( Nijconj ) ) < TOL ) )
      return Tf[(int) Nij][(int) Nji][(int) Nijconj];

    for( i = 0; i < Tijdom[0][1]; i++ )
      if ( Nij >= i && Nij <= i + 1 )
        x = i;
    for( i = 0; i < Tijdom[1][1]; i++ )
      if ( Nji >= i && Nji <= i + 1 )
        y = i;
    for( i = 0; i < Tijdom[2][1]; i++ )
      if ( Nijconj >= i && Nijconj <= i + 1 )
        z = i;
    for( i = 0; i < 64; i++ )
      coeffs[i] = Tijc[x][y][z][i];

    return Sptricubic( Nij, Nji, Nijconj, coeffs );
  }

  double AIREBO::Sp5th( double x, double *coeffs ) {
    double f;
    const double x2 = x * x;
    const double x3 = x2 * x;

    f = coeffs[0];
    f += coeffs[1] * x;
    f += coeffs[2] * x2;
    f += coeffs[3] * x3;
    f += coeffs[4] * x2 * x2;
    f += coeffs[5] * x2 * x3;

    return f;
  }

  double AIREBO::Spbicubic( double x, double y, double *coeffs ) {
    double f = 0.0, xn; // yn values can actually be cached

    double yn[DIM];
    yn[0] = 1.0;
    for( int i = 1; i < DIM; i++ )
      yn[i] = yn[i - 1] * y;

    xn = 1.0;
    for( int i_x = 0; i_x < DIM; i_x++ ) {
      for( int i_y = 0; i_y < DIM; i_y++ )
        f += coeffs[i_x * DIM + i_y] * xn * yn[i_y];
      xn *= x;
    }

    return f;
  }

  double AIREBO::Sptricubic( double x, double y, double z, double *coeffs ) {
    double f = 0.0, xn; // yn & zn values can actually be cached

    double yn[DIM], zn[DIM];
    yn[0] = 1.0;
    zn[0] = 1.0;
    for( int i = 1; i < DIM; i++ ) {
      yn[i] = yn[i - 1] * y;
      zn[i] = zn[i - 1] * z;
    }

    xn = 1.0;
    for( int i_x = 0; i_x < DIM; i_x++ ) {
      for( int i_y = 0; i_y < DIM; i_y++ )
        for( int i_z = 0; i_z < DIM; i_z++ )
          f += coeffs[DIM * DIM * i_x + DIM * i_y + i_z] * xn * yn[i_y] * zn[i_z];
      xn *= x;
    }

    return f;
  }

}