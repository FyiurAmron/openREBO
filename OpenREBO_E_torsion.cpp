/*
    File:   OpenREBO_E_Torsion.cpp
    Author: vaxquis
    based on LAMMPS implementation of AIREBO potential
 */

#include "OpenREBO.h"

#ifndef LEAN_REBO

namespace OpenREBO {

  double ForceField::E_Torsion( ) {
    double eTorsion = 0.0;

    int itype;
    int j, jtype;
    int k, ktype;
    int l, ltype;

    int REBO_neighbours_num_i;
    int *REBO_neighbours_list_i;
    vec3d *REBO_neighbours_bonds_i;

    int REBO_neighbours_num_j;
    int *REBO_neighbours_list_j;
    vec3d *REBO_neighbours_bonds_j;

    double rji_vec[3], rji, rji_sq;
    double rij_vec[3], rij, rij_sq;
    double rki_vec[3], rki, rki_sq;
    double rlj_vec[3], rlj, rlj_sq;
    double rkj_vec[3], rkj_sq;
    double rli_vec[3], rli_sq;
    double wji, wki, wlj;

    double cos321, sin321;
    double cos234, sin234;
    double cross321[3], cross234[3];

    double cwnum, cwnom, cw, cw2;

    double tspjik, tspijl;

    double Vtors;

    for( int i = 0; i < natoms; i++ ) {
      itype = type[i];

      if ( itype != 0 )
        continue;

      REBO_neighbours_num_i = REBO_neighbours_num[i];
      REBO_neighbours_list_i = REBO_neighbours_list[i];
      REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];

      for( int jj = 0; jj < REBO_neighbours_num_i; jj++ ) {
        j = REBO_neighbours_list_i[jj];

        if ( i < j )
          continue;

        jtype = type[j];

        if ( jtype != 0 )
          continue;

        rji_vec[0] = REBO_neighbours_bonds_i[jj].x;
        rji_vec[1] = REBO_neighbours_bonds_i[jj].y;
        rji_vec[2] = REBO_neighbours_bonds_i[jj].z;
        rji = REBO_neighbours_bonds_i[jj].r;
        rji_sq = REBO_neighbours_bonds_i[jj].r_sq;

        neg( rij_vec, rji_vec );
        rij = rji;
        rij_sq = rji_sq;

        wji = SpRC( rji, itype, jtype );

        for( int kk = 0; kk < REBO_neighbours_num_i; kk++ ) {
          k = REBO_neighbours_list_i[kk];

          ktype = type[k];

          if ( k == j )
            continue;

          rki_vec[0] = REBO_neighbours_bonds_i[kk].x;
          rki_vec[1] = REBO_neighbours_bonds_i[kk].y;
          rki_vec[2] = REBO_neighbours_bonds_i[kk].z;
          rki = REBO_neighbours_bonds_i[kk].r;
          rki_sq = REBO_neighbours_bonds_i[kk].r_sq;

          cos321 = -cos_theta_clamp( rki_vec, rij_vec, rki, rij );
          sin321 = sqrt( 1.0 - cos321 * cos321 );

          if ( sin321 < TOL )
            continue;

          diff( rkj_vec, rki_vec, rji_vec );
          rkj_sq = length_sq( rkj_vec );

          wki = SpRC( rki, itype, ktype );

          tspjik = Sp2th( 0.5 * ( rij_sq + rki_sq - rkj_sq ) / ( rij * rki ) );

          REBO_neighbours_num_j = REBO_neighbours_num[j];
          REBO_neighbours_list_j = REBO_neighbours_list[j];
          REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];

          for( int ll = 0; ll < REBO_neighbours_num_j; ll++ ) {
            l = REBO_neighbours_list_j[ll];

            ltype = type[l];

            if ( l == i || l == k )
              continue;

            rlj_vec[0] = REBO_neighbours_bonds_j[ll].x;
            rlj_vec[1] = REBO_neighbours_bonds_j[ll].y;
            rlj_vec[2] = REBO_neighbours_bonds_j[ll].z;
            rlj = REBO_neighbours_bonds_j[ll].r;
            rlj_sq = REBO_neighbours_bonds_j[ll].r_sq;

            cos234 = cos_theta_clamp( rij_vec, rlj_vec, rij, rlj );
            sin234 = sqrt( 1.0 - cos234 * cos234 );

            if ( sin234 < TOL )
              continue;

            wlj = SpRC( rlj, jtype, ltype );

            sum( rli_vec, rji_vec, rlj_vec );
            rli_sq = length_sq( rli_vec );

            tspijl = Sp2th( 0.5 * ( rji_sq + rlj_sq - rli_sq ) / ( rji * rlj ) );

            cross( cross321, rij_vec, rki_vec );
            cross( cross234, rji_vec, rlj_vec );

            cwnum = dot( cross321, cross234 );
            cwnom = rki * rlj * rij * rij * sin321 * sin234;

            cw = 0.5 - 0.5 * cwnum / cwnom;
            cw2 = cw * cw;
            cw2 *= cw2;
            cw2 *= cw;
            Vtors = epsilonT[ktype][ltype] * ( 256.0 / 405.0 * cw2 - 0.1 );

            eTorsion += Vtors * wki * wji * wlj * ( 1.0 - tspjik ) * ( 1.0 - tspijl );
          }
        }
      }
    }
    return eTorsion;
  }

}

#endif