// based on LAMMPS implementation of AIREBO force field

#include "airebo_force_field.h"

namespace AIREBO {

  void ForceField::E_TORSION( ) {
    int itype;
    int j, jtype;
    int k, ktype;
    int l, ltype;
    int jj, kk, ll;

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
    double costmp;
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

      for( jj = 0; jj < REBO_neighbours_num_i; jj++ ) {
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

        rij_vec[0] = -rji_vec[0];
        rij_vec[1] = -rji_vec[1];
        rij_vec[2] = -rji_vec[2];
        rij = rji;
        rij_sq = rji_sq;

        wji = Sp( rji, rcMin[itype][jtype], rcMax[itype][jtype] );

        for( kk = 0; kk < REBO_neighbours_num_i; kk++ ) {
          k = REBO_neighbours_list_i[kk];

          ktype = type[k];

          if ( k == j )
            continue;

          rki_vec[0] = REBO_neighbours_bonds_i[kk].x;
          rki_vec[1] = REBO_neighbours_bonds_i[kk].y;
          rki_vec[2] = REBO_neighbours_bonds_i[kk].z;
          rki = REBO_neighbours_bonds_i[kk].r;
          rki_sq = REBO_neighbours_bonds_i[kk].r_sq;

          cos321 = -( ( rki_vec[0] * rij_vec[0] ) + ( rki_vec[1] * rij_vec[1] ) + ( rki_vec[2] * rij_vec[2] ) ) / ( rki * rij );
          cos321 = min( cos321, 1.0 );
          cos321 = max( cos321, -1.0 );
          sin321 = sqrt( 1.0 - cos321 * cos321 );

          if ( sin321 < TOL )
            continue;

          rkj_vec[0] = rki_vec[0] - rji_vec[0];
          rkj_vec[1] = rki_vec[1] - rji_vec[1];
          rkj_vec[2] = rki_vec[2] - rji_vec[2];

          rkj_sq = rkj_vec[0] * rkj_vec[0] + rkj_vec[1] * rkj_vec[1] + rkj_vec[2] * rkj_vec[2];

          wki = Sp( rki, rcMin[itype][ktype], rcMax[itype][ktype] );

          costmp = 0.5 * ( rij_sq + rki_sq - rkj_sq ) / rij / rki;
          tspjik = Sp2th( costmp );

          REBO_neighbours_num_j = REBO_neighbours_num[j];
          REBO_neighbours_list_j = REBO_neighbours_list[j];
          REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];

          for( ll = 0; ll < REBO_neighbours_num_j; ll++ ) {
            l = REBO_neighbours_list_j[ll];

            ltype = type[l];

            if ( ( l == i ) || ( l == k ) )
              continue;

            rlj_vec[0] = REBO_neighbours_bonds_j[ll].x;
            rlj_vec[1] = REBO_neighbours_bonds_j[ll].y;
            rlj_vec[2] = REBO_neighbours_bonds_j[ll].z;
            rlj = REBO_neighbours_bonds_j[ll].r;
            rlj_sq = REBO_neighbours_bonds_j[ll].r_sq;

            cos234 = ( rij_vec[0] * rlj_vec[0] + rij_vec[1] * rlj_vec[1] + rij_vec[2] * rlj_vec[2] ) / ( rij * rlj );
            cos234 = min( cos234, 1.0 );
            cos234 = max( cos234, -1.0 );
            sin234 = sqrt( 1.0 - cos234 * cos234 );

            if ( sin234 < TOL )
              continue;

            wlj = Sp( rlj, rcMin[jtype][ltype], rcMax[jtype][ltype] );

            rli_vec[0] = rji_vec[0] + rlj_vec[0];
            rli_vec[1] = rji_vec[1] + rlj_vec[1];
            rli_vec[2] = rji_vec[2] + rlj_vec[2];
            rli_sq = rli_vec[0] * rli_vec[0] + rli_vec[1] * rli_vec[1] + rli_vec[2] * rli_vec[2];

            costmp = 0.5 * ( rji_sq + rlj_sq - rli_sq ) / rji / rlj;
            tspijl = Sp2th( costmp );

            cross321[0] = ( rij_vec[1] * rki_vec[2] ) - ( rij_vec[2] * rki_vec[1] );
            cross321[1] = ( rij_vec[2] * rki_vec[0] ) - ( rij_vec[0] * rki_vec[2] );
            cross321[2] = ( rij_vec[0] * rki_vec[1] ) - ( rij_vec[1] * rki_vec[0] );

            cross234[0] = ( rji_vec[1] * rlj_vec[2] ) - ( rji_vec[2] * rlj_vec[1] );
            cross234[1] = ( rji_vec[2] * rlj_vec[0] ) - ( rji_vec[0] * rlj_vec[2] );
            cross234[2] = ( rji_vec[0] * rlj_vec[1] ) - ( rji_vec[1] * rlj_vec[0] );

            cwnum = ( cross321[0] * cross234[0] ) + ( cross321[1] * cross234[1] ) + ( cross321[2] * cross234[2] );
            cwnom = rki * rlj * rij * rij * sin321 * sin234;

            cw = 0.5 - 0.5 * cwnum / cwnom;
            cw2 = cw * cw;
            cw2 *= cw2;
            cw2 *= cw;
            Vtors = epsilonT[ktype][ltype] * ( 256.0 / 405.0 * cw2 - 0.1 );

            energy_torsion += Vtors * wki * wji * wlj * ( 1.0 - tspjik ) * ( 1.0 - tspijl );
          }
        }
      }
    }
  }

}