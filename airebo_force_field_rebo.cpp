// written by Szymon Winczewski
// based on LAMMPS implementation of AIREBO force field

#include "airebo_force_field.h"
#include <cstring>

using namespace std;

void AIREBOForceField::REBO_neighbours( ) {
  int i, itype;
  int j, jtype;
  int jj;

  int REBO_neighbours_num_i;
  int *REBO_neighbours_list_i;
  vec3d *REBO_neighbours_bonds_i;

  int neighbours_num_i;
  int *neighbours_list_i;
  vec3d *neighbours_bonds_i;

  double tmp_r_sq;

  for( i = 0; i < natoms; i++ ) {
    REBO_neighbours_num_i = 0;
    REBO_neighbours_list_i = REBO_neighbours_list[i];
    REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];
    itype = type[i];
    //nC[i] = nH[i] = 0.0; // already set by nC/nH new[]()

    neighbours_num_i = neighbours_num[i];
    neighbours_list_i = neighbours_list[i];
    neighbours_bonds_i = neighbours_bonds[i];

    for( jj = 0; jj < neighbours_num_i; jj++ ) {
      j = neighbours_list_i[jj];
      jtype = type[j];

      tmp_r_sq = neighbours_bonds_i[jj].r_sq;

      if ( tmp_r_sq >= rcMaxSq[itype][jtype] )
        continue;

      REBO_neighbours_list_i[REBO_neighbours_num_i] = j;
      memcpy( REBO_neighbours_bonds_i[REBO_neighbours_num_i], neighbours_bonds_i[jj], sizeof (vec3d ) )
      REBO_neighbours_bonds_i[REBO_neighbours_num_i].x = neighbours_bonds_i[jj].x;
      REBO_neighbours_bonds_i[REBO_neighbours_num_i].y = neighbours_bonds_i[jj].y;
      REBO_neighbours_bonds_i[REBO_neighbours_num_i].z = neighbours_bonds_i[jj].z;
      REBO_neighbours_bonds_i[REBO_neighbours_num_i].r = neighbours_bonds_i[jj].r;
      REBO_neighbours_bonds_i[REBO_neighbours_num_i].r_sq = tmp_r_sq;
      REBO_neighbours_num_i++;

      if ( jtype == 0 )
        //nC[i] += Sp( sqrt( tmp_r_sq ), rcmin[itype][jtype], rcmax[itype][jtype] );
        nC[i] += SpRC( sqrt( tmp_r_sq ), itype, jtype );
      else
        nH[i] += SpRC( sqrt( tmp_r_sq ), itype, jtype );
      //nH[i] += Sp( sqrt( tmp_r_sq ), rcmin[itype][jtype], rcmax[itype][jtype] );
    }
    // TODO delegate check to nlists loader
    if ( REBO_neighbours_num_i > max_number_of_REBO_neighbours ) {
      cout << "error in AIREBOForceField::REBO_neighbours(): maximum number of REBO neighbours reached!" << endl;
      exit( 0 );
    }

    REBO_neighbours_num[i] = REBO_neighbours_num_i;
  }
}

void AIREBOForceField::E_REBO( ) {
  int i, j;
  int jj;
  int itype, jtype;
  int REBO_neighbours_num_i;
  int *REBO_neighbours_list_i;
  vec3d *REBO_neighbours_bonds_i;
  double rji, rji_vec[3];
  double wij;
  double Qij, Aij, alphaij;
  double VR, VA;
  double term;
  int m;
  double bij;

  for( i = 0; i < natoms; i++ ) {
    itype = type[i];
    REBO_neighbours_num_i = REBO_neighbours_num[i];
    REBO_neighbours_list_i = REBO_neighbours_list[i];
    REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];

    for( jj = 0; jj < REBO_neighbours_num_i; jj++ ) {
      j = REBO_neighbours_list_i[jj];

      if ( i > j )
        continue;

      jtype = type[j];

      rji = REBO_neighbours_bonds_i[jj].r;

      wij = Sp( rji, rcMin[itype][jtype], rcMax[itype][jtype] );
      if ( wij <= TOL )
        continue;

      rji_vec[0] = REBO_neighbours_bonds_i[jj].x;
      rji_vec[1] = REBO_neighbours_bonds_i[jj].y;
      rji_vec[2] = REBO_neighbours_bonds_i[jj].z;

      Qij = Q[itype][jtype];
      Aij = A[itype][jtype];
      alphaij = alpha[itype][jtype];

      VR = wij * ( 1.0 + ( Qij / rji ) ) * Aij * exp( -alphaij * rji );

      VA = 0.0;
      for( m = 0; m < 3; m++ ) {
        term = -wij * BIJc[itype][jtype][m] * exp( -Beta[itype][jtype][m] * rji );
        VA += term;
      }

      bij = bondorder( i, j, rji_vec, rji );
      energy_rebo += ( VR + bij * VA );
    }
  }
}

double AIREBOForceField::bondorder( int i, int j, double rji_vec[3], double rji ) {
  int atomi, atomj;
  int itype, jtype;
  int k, atomk, ktype;
  int l, atoml, ltype;

  int REBO_neighbours_num_i;
  int *REBO_neighbours_list_i;
  vec3d *REBO_neighbours_bonds_i;
  int REBO_neighbours_num_j;
  int *REBO_neighbours_list_j;
  vec3d *REBO_neighbours_bonds_j;

  double rij_vec[3], rij, rij_sq;
  double rki_vec[3], rki, rki_sq;
  double rlj_vec[3], rlj, rlj_sq;
  double rkj_vec[3], rkj_sq;
  double rli_vec[3], rli_sq;

  double lamdajik, lamdaijl;
  double wji, wki, wlj;
  double cosjik, cosijl;
  double Nki, Nlj;

  double g;
  double PijS, PjiS;
  double bij;

  double NijC, NijH, NjiC, NjiH;

  double NconjtmpI, NconjtmpJ, Nijconj;

  double pij, pji;
  double piRC, Tij, Etmp;

  double cos321, sin321;
  double cos234, sin234;

  double costmp;
  double tspjik, tspijl;
  double cross321[3], cross234[3];

  double cwnum, cwnom;
  double om1234;

  atomi = i;
  atomj = j;
  itype = type[i];
  jtype = type[j];

  wji = Sp( rji, rcMin[itype][jtype], rcMax[itype][jtype] );

  NijC = nC[i] - ( wji * kronecker( jtype, 0 ) );
  NijH = nH[i] - ( wji * kronecker( jtype, 1 ) );
  NjiC = nC[j] - ( wji * kronecker( itype, 0 ) );
  NjiH = nH[j] - ( wji * kronecker( itype, 1 ) );

  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;

  REBO_neighbours_num_i = REBO_neighbours_num[i];
  REBO_neighbours_list_i = REBO_neighbours_list[i];
  REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];
  for( k = 0; k < REBO_neighbours_num_i; k++ ) {
    atomk = REBO_neighbours_list_i[k];
    if ( atomk != atomj ) {
      ktype = type[atomk];
      rki_vec[0] = REBO_neighbours_bonds_i[k].x;
      rki_vec[1] = REBO_neighbours_bonds_i[k].y;
      rki_vec[2] = REBO_neighbours_bonds_i[k].z;
      rki = REBO_neighbours_bonds_i[k].r;

      lamdajik = 4.0 * kronecker( itype, 1 ) * ( ( rho[ktype][1] - rki ) - ( rho[jtype][1] - rji ) );
      wki = Sp( rki, rcMin[itype][ktype], rcMax[itype][ktype] );
      Nki = nC[atomk] - ( wki * kronecker( itype, 0 ) ) + nH[atomk] - ( wki * kronecker( itype, 1 ) );
      cosjik = ( ( rji_vec[0] * rki_vec[0] ) + ( rji_vec[1] * rki_vec[1] ) + ( rji_vec[2] * rki_vec[2] ) ) / ( rji * rki );
      cosjik = min( cosjik, 1.0 );
      cosjik = max( cosjik, -1.0 );
      g = gSpline( cosjik, ( NijC + NijH ), itype );
      Etmp = Etmp + ( wki * g * exp( lamdajik ) );
      NconjtmpI = NconjtmpI + ( kronecker( ktype, 0 ) * wki * Sp( Nki, Nmin, Nmax ) );
    }
  }

  PijS = 0.0;
  PijS = PijSpline( NijC, NijH, itype, jtype );
  pij = pow( 1.0 + Etmp + PijS, -0.5 );

  Etmp = 0.0;

  REBO_neighbours_num_j = REBO_neighbours_num[j];
  REBO_neighbours_list_j = REBO_neighbours_list[j];
  REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];
  for( l = 0; l < REBO_neighbours_num_j; l++ ) {
    atoml = REBO_neighbours_list_j[l];
    if ( atoml != atomi ) {
      ltype = type[atoml];
      rlj_vec[0] = REBO_neighbours_bonds_j[l].x;
      rlj_vec[1] = REBO_neighbours_bonds_j[l].y;
      rlj_vec[2] = REBO_neighbours_bonds_j[l].z;
      rlj = REBO_neighbours_bonds_j[l].r;

      lamdaijl = 4.0 * kronecker( jtype, 1 ) * ( ( rho[ltype][1] - rlj ) - ( rho[itype][1] - rji ) );
      wlj = Sp( rlj, rcMin[jtype][ltype], rcMax[jtype][ltype] );
      Nlj = nC[atoml] - ( wlj * kronecker( jtype, 0 ) ) + nH[atoml] - ( wlj * kronecker( jtype, 1 ) );
      cosijl = -1.0 * ( ( rji_vec[0] * rlj_vec[0] ) + ( rji_vec[1] * rlj_vec[1] ) + ( rji_vec[2] * rlj_vec[2] ) ) / ( rji * rlj );
      cosijl = min( cosijl, 1.0 );
      cosijl = max( cosijl, -1.0 );

      g = gSpline( cosijl, NjiC + NjiH, jtype );
      Etmp = Etmp + ( wlj * g * exp( lamdaijl ) );
      NconjtmpJ = NconjtmpJ + ( kronecker( ltype, 0 ) * wlj * Sp( Nlj, Nmin, Nmax ) );
    }
  }

  PjiS = 0.0;
  PjiS = PijSpline( NjiC, NjiH, jtype, itype );
  pji = pow( 1.0 + Etmp + PjiS, -0.5 );

  Nijconj = 1.0 + ( NconjtmpI * NconjtmpI ) + ( NconjtmpJ * NconjtmpJ );
  piRC = piRCSpline( NijC + NijH, NjiC + NjiH, Nijconj, itype, jtype );

  Tij = 0.0;
  if ( ( itype == 0 ) && ( jtype == 0 ) )
    Tij = TijSpline( ( NijC + NijH ), ( NjiC + NjiH ), Nijconj );
  Etmp = 0.0;

  if ( fabs( Tij ) > TOL ) {
    rij_vec[0] = -rji_vec[0];
    rij_vec[1] = -rji_vec[1];
    rij_vec[2] = -rji_vec[2];
    rij = rji;
    rij_sq = rij * rij;

    REBO_neighbours_num_i = REBO_neighbours_num[i];
    REBO_neighbours_list_i = REBO_neighbours_list[i];
    REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];
    for( k = 0; k < REBO_neighbours_num_i; k++ ) {
      atomk = REBO_neighbours_list_i[k];
      ktype = type[atomk];
      if ( atomk != atomj ) {
        rki_vec[0] = REBO_neighbours_bonds_i[k].x;
        rki_vec[1] = REBO_neighbours_bonds_i[k].y;
        rki_vec[2] = REBO_neighbours_bonds_i[k].z;
        rki = REBO_neighbours_bonds_i[k].r;
        rki_sq = REBO_neighbours_bonds_i[k].r_sq;

        cos321 = -1.0 * ( ( rki_vec[0] * rij_vec[0] ) + ( rki_vec[1] * rij_vec[1] ) + ( rki_vec[2] * rij_vec[2] ) ) / ( rki * rij );
        cos321 = min( cos321, 1.0 );
        cos321 = max( cos321, -1.0 );

        sin321 = sqrt( 1.0 - cos321 * cos321 );
        if ( sin321 != 0.0 ) {
          rkj_vec[0] = rki_vec[0] - rji_vec[0];
          rkj_vec[1] = rki_vec[1] - rji_vec[1];
          rkj_vec[2] = rki_vec[2] - rji_vec[2];
          rkj_sq = ( rkj_vec[0] * rkj_vec[0] ) + ( rkj_vec[1] * rkj_vec[1] ) + ( rkj_vec[2] * rkj_vec[2] );
          wki = Sp( rki, rcMin[itype][ktype], rcMaxP[itype][ktype] );
          costmp = 0.5 * ( rij_sq + rki_sq - rkj_sq ) / rij / rki;
          tspjik = Sp2( costmp, thmin, thmax );

          REBO_neighbours_num_j = REBO_neighbours_num[j];
          REBO_neighbours_list_j = REBO_neighbours_list[j];
          REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];
          for( l = 0; l < REBO_neighbours_num_j; l++ ) {
            atoml = REBO_neighbours_list_j[l];
            ltype = type[atoml];
            if ( !( ( atoml == atomi ) || ( atoml == atomk ) ) ) {
              rlj_vec[0] = REBO_neighbours_bonds_j[l].x;
              rlj_vec[1] = REBO_neighbours_bonds_j[l].y;
              rlj_vec[2] = REBO_neighbours_bonds_j[l].z;
              rlj = REBO_neighbours_bonds_j[l].r;
              rlj_sq = REBO_neighbours_bonds_j[l].r_sq;

              cos234 = ( rij_vec[0] * rlj_vec[0] + rij_vec[1] * rlj_vec[1] + rij_vec[2] * rlj_vec[2] ) / ( rij * rlj );
              cos234 = min( cos234, 1.0 );
              cos234 = max( cos234, -1.0 );
              sin234 = sqrt( 1.0 - cos234 * cos234 );

              if ( sin234 != 0.0 ) {
                wlj = Sp( rlj, rcMin[jtype][ltype], rcMaxP[jtype][ltype] );
                rli_vec[0] = rji_vec[0] + rlj_vec[0];
                rli_vec[1] = rji_vec[1] + rlj_vec[1];
                rli_vec[2] = rji_vec[2] + rlj_vec[2];
                rli_sq = ( rli_vec[0] * rli_vec[0] ) + ( rli_vec[1] * rli_vec[1] ) + ( rli_vec[2] * rli_vec[2] );

                costmp = 0.5 * ( rij_sq + rlj_sq - rli_sq ) / rij / rlj;
                tspijl = Sp2( costmp, thmin, thmax );

                cross321[0] = ( rij_vec[1] * rki_vec[2] ) - ( rij_vec[2] * rki_vec[1] );
                cross321[1] = ( rij_vec[2] * rki_vec[0] ) - ( rij_vec[0] * rki_vec[2] );
                cross321[2] = ( rij_vec[0] * rki_vec[1] ) - ( rij_vec[1] * rki_vec[0] );
                cross234[0] = ( rji_vec[1] * rlj_vec[2] ) - ( rji_vec[2] * rlj_vec[1] );
                cross234[1] = ( rji_vec[2] * rlj_vec[0] ) - ( rji_vec[0] * rlj_vec[2] );
                cross234[2] = ( rji_vec[0] * rlj_vec[1] ) - ( rji_vec[1] * rlj_vec[0] );

                cwnum = ( cross321[0] * cross234[0] ) + ( cross321[1] * cross234[1] ) + ( cross321[2] * cross234[2] );
                cwnom = rki * rlj * rji * rji * sin321 * sin234;
                om1234 = cwnum / cwnom;

                Etmp += ( ( 1.0 - square( om1234 ) ) * wki * wlj ) * ( 1.0 - tspjik ) * ( 1.0 - tspijl );
              }
            }
          }
        }
      }
    }
  }

  bij = ( 0.5 * ( pij + pji ) ) + piRC + ( Tij * Etmp );
  return bij;
}
