/*
    File:   splines.cpp
    Author: vaxquis
    based on LAMMPS implementation of AIREBO potential
 */

#include "OpenREBO.h"

namespace OpenREBO {
  double AIREBO::Sptricubic( double x, double y, double z, double coeffs[SPL_DIM*SPL_DIM*SPL_DIM] ) {
    double f = 0.0;

    for( int i_x = SPL_DIM - 1; i_x > 0; i_x-- ) {
      for( int i_y = SPL_DIM - 1; i_y > 0; i_y-- ) {
        for( int i_z = SPL_DIM - 1; i_z > 0; i_z-- ) {
          f += coeffs[SPL_DIM * SPL_DIM * i_x + SPL_DIM * i_y + i_z];
          f *= z;
        }
        f += coeffs[SPL_DIM * SPL_DIM * i_x + SPL_DIM * i_y];
        f *= y;
      }
      for( int i_z = SPL_DIM - 1; i_z > 0; i_z-- ) {
        f += coeffs[SPL_DIM * SPL_DIM * i_x + i_z];
        f *= z;
      }
      f += coeffs[SPL_DIM * SPL_DIM * i_x];
      f *= x;
    }

    for( int i_y = SPL_DIM - 1; i_y > 0; i_y-- ) {
      for( int i_z = SPL_DIM - 1; i_z > 0; i_z-- ) {
        f += coeffs[SPL_DIM * i_y + i_z];
        f *= z;
      }
      f += coeffs[SPL_DIM * i_y];
      f *= y;
    }
    for( int i_z = SPL_DIM - 1; i_z > 0; i_z-- ) {
      f += coeffs[i_z];
      f *= z;
    }

    return f + coeffs[0];
  }

  double AIREBO::Spbicubic( double x, double y, double coeffs[SPL_DIM*SPL_DIM] ) {
    double f = 0.0;

    for( int i_x = SPL_DIM - 1; i_x > 0; i_x-- ) {
      for( int i_y = SPL_DIM - 1; i_y > 0; i_y-- ) {
        f += coeffs[SPL_DIM * i_x + i_y];
        f *= y;
      }
      f += coeffs[SPL_DIM * i_x];
      f *= x;
    }
    for( int i_z = SPL_DIM - 1; i_z > 0; i_z-- ) {
      f += coeffs[i_z];
      f *= y;
    }

    return f + coeffs[0];
  }

  double AIREBO::Sp5th( double x, double coeffs[6] ) {
    double f = 0.0;
    for( int i = 5; i > 0; i-- ) {
      f += coeffs[i];
      f *= x;
    }
    return f + coeffs[0];
  }

  double AIREBO::gSpline( double costh, double Nij, int typei ) {
    return ( typei == 0 ) ? gSplineC( costh, Nij ) : gSplineH( costh );
  }

  double AIREBO::gSplineC( double costh, double Nij ) {
    if ( Nij <= NCmin ) {
      if ( costh <= gCdom[0] )
        return gSplineC1_low;
      if ( costh >= gCdom[4] )
        return gSplineC1_hi;

      for( int i = 1; i < 4; i++ )
        if ( costh <= gCdom[i] )
          return Sp5th( costh, gC1[i - 1] );
      return Sp5th( costh, gC1[3] );
    }

    if ( Nij >= NCmax ) {
      if ( costh <= gCdom[0] )
        return gSplineC2_low;
      if ( costh >= gCdom[4] )
        return gSplineC2_hi;

      for( int i = 1; i < 4; i++ )
        if ( costh <= gCdom[i] )
          return Sp5th( costh, gC2[i - 1] );
      return Sp5th( costh, gC2[3] );
    }

    // in both bounds
    double g2;

    if ( costh <= gCdom[0] )
      return gSplineC2_low + SpNC( Nij ) * gSplineC_low_delta;
    if ( costh >= gCdom[4] )
      return gSplineC2_hi + SpNC( Nij ) * gSplineC_hi_delta;

    for( int i = 1; i < 4; i++ )
      if ( costh <= gCdom[i] ) {
        i--;
        g2 = Sp5th( costh, gC2[i] );
        return g2 + SpNC( Nij ) * ( Sp5th( costh, gC1[i] ) - g2 );
      }
    g2 = Sp5th( costh, gC2[3] );
    return g2 + SpNC( Nij ) * ( Sp5th( costh, gC1[3] ) - g2 );
  }

  double AIREBO::gSplineH( double costh ) { // central atom is Hydrogen
    if ( costh < gHdom[0] )
      return gSplineH_low;
    if ( costh > gHdom[3] )
      return gSplineH_hi;

    for( int i = 1; i < 3; i++ )
      if ( costh <= gHdom[i] )
        return Sp5th( costh, gH[i - 1] );
    return Sp5th( costh, gH[2] );
  }

  double AIREBO::PijSplineC_only( double NijC ) {
    if ( NijC < pCCdom[0][0] )
      NijC = pCCdom[0][0];
    else if ( NijC > pCCdom[0][1] )
      NijC = pCCdom[0][1];

    int iNijC = (int) NijC;
    return ( ( NijC - iNijC < TOL ) ?
            PCCf[iNijC][0]
            : pCC[iNijC][0][0] );
  }

  double AIREBO::PijSpline( double NijC, double NijH, int typei, int typej ) {
    if ( typei == 1 )
      return 0.0;

    // if the inputs are out of bounds set them back to a point in bounds
    // note: lower bound is 0 (by pCXdom values),
    // otherwise the code below gets broken by out-of-bounds array access
    if ( typej == 0 ) {
      if ( NijC < pCCdom[0][0] )
        NijC = pCCdom[0][0];
      else if ( NijC > pCCdom[0][1] )
        NijC = pCCdom[0][1];
      if ( NijH < pCCdom[1][0] )
        NijH = pCCdom[1][0];
      else if ( NijH > pCCdom[1][1] )
        NijH = pCCdom[1][1];

      return ( ( NijC - (int) NijC < TOL ) && ( NijH - int( NijH ) < TOL ) ) ?
              PCCf[(int) NijC][(int) NijH]
              : Spbicubic( NijC, NijH, pCC[(int) NijC ][(int) NijH ] );
    }

    // else ( typej == 1 )
    if ( NijC < pCHdom[0][0] )
      NijC = pCHdom[0][0];
    else if ( NijC > pCHdom[0][1] )
      NijC = pCHdom[0][1];
    if ( NijH < pCHdom[1][0] )
      NijH = pCHdom[1][0];
    else if ( NijH > pCHdom[1][1] )
      NijH = pCHdom[1][1];

    return ( ( NijC - (int) NijC < TOL ) && ( NijH - int( NijH ) < TOL ) ) ?
            PCHf[(int) NijC][(int) NijH]
            : Spbicubic( NijC, NijH, pCH[(int) NijC ][(int) NijH ] );
  }

  double AIREBO::TricubicSplineMain( double Nij, double Nji, double Nijconj,
          double dom[3][2], double coeff_f[5][5][11], double coeff[SPL_DIM][SPL_DIM][9][SPL_DIM*SPL_DIM*SPL_DIM] ) {
    if ( Nij < dom[0][0] )
      Nij = dom[0][0];
    else if ( Nij > dom[0][1] )
      Nij = dom[0][1];
    if ( Nji < dom[1][0] )
      Nji = dom[1][0];
    else if ( Nji > dom[1][1] )
      Nji = dom[1][1];
    if ( Nijconj < dom[2][0] )
      Nijconj = dom[2][0];
    else if ( Nijconj > dom[2][1] )
      Nijconj = dom[2][1];

    int iNij = (int) Nij, iNji = (int) Nji, iNijconj = (int) Nijconj;
    if ( ( Nij - iNij < TOL ) && ( Nji - iNji < TOL ) && ( Nijconj - iNijconj < TOL ) )
      return coeff_f[iNij][iNji][iNijconj];

    // NOTE: ain't round() better here?
    return Sptricubic( Nij, Nji, Nijconj, coeff[iNij][iNji][iNijconj] );
  }

  double AIREBO::piRCSplineCC( double Nij, double Nji, double Nijconj ) {
    return TricubicSplineMain( Nij, Nji, Nijconj, piCCdom, piCCf, piCC );
  }

  double AIREBO::piRCSplineCH( double Nij, double Nji, double Nijconj ) {
    return TricubicSplineMain( Nij, Nji, Nijconj, piCHdom, piCHf, piCH );
  }

  double AIREBO::piRCSplineHH( double Nij, double Nji, double Nijconj ) {
    return TricubicSplineMain( Nij, Nji, Nijconj, piHHdom, piHHf, piHH );
  }

  double AIREBO::piRCSpline( double Nij, double Nji, double Nijconj, int typei, int typej ) {
    return (typei == 0 ) ?
            ( ( typej == 0 ) ? piRCSplineCC( Nij, Nji, Nijconj ) : piRCSplineCH( Nij, Nji, Nijconj ) )
            : ( ( typej == 0 ) ? piRCSplineCH( Nij, Nji, Nijconj ) : piRCSplineHH( Nij, Nji, Nijconj ) );
  }

  double AIREBO::TijSpline( double Nij, double Nji, double Nijconj ) {
    return TricubicSplineMain( Nij, Nji, Nijconj, Tijdom, Tf, Tijc );
  }
}