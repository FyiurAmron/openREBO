/*
    File:   E_LJ.cpp
    Author: vaxquis
    based on LAMMPS implementation of AIREBO potential
 */

#include "OpenREBO.h"

#ifndef LEAN_REBO

namespace OpenREBO {

  double ForceField::E_LJ( ) {
    double eLJ = 0.0;
    double rljmin, rljmax;
    double sigcut, sigmin, sigwid;
    int i, j, k, m;
    int jj, kk, mm;
    int itype, jtype, ktype, mtype;

    int neighbours_num_i;
    int *neighbours_list_i;
    vec3d *neighbours_bonds_i;

    int REBO_neighbours_num_i;
    int *REBO_neighbours_list_i;
    vec3d *REBO_neighbours_bonds_i;

    int REBO_neighbours_num_k;
    int *REBO_neighbours_list_k;
    vec3d *REBO_neighbours_bonds_k;

    double best;
    bool testpath, done;

    double rji, rji_sq, rji_vec[3];
    double rki, rki_sq, rki_vec[3];
    double rkj, rkj_sq, rkj_vec[3];
    double rmk, rmk_sq, rmk_vec[3];
    double rmj, rmj_sq, rmj_vec[3];

    double wki, wkj, wmk, wmj;

    double cij;
    double slw;
    double drji;
    double swidth, tee, tee2;
    double r2inv, r6inv;
    double vdw, VLJ, Str, VA, Stb;
    double delscale[3];

    rljmin = 0.0;
    rljmax = 0.0;
    sigcut = 0.0;
    sigmin = 0.0;
    sigwid = 0.0;

    for( i = 0; i < natoms; i++ ) {
      itype = type[i];

      neighbours_num_i = neighbours_num[i];
      neighbours_list_i = neighbours_list[i];
      neighbours_bonds_i = neighbours_bonds[i];

      for( jj = 0; jj < neighbours_num_i; jj++ ) {
        j = neighbours_list_i[jj];

        if ( i < j )
          continue;

        jtype = type[j];

        rji_sq = neighbours_bonds_i[jj].r_sq;
        if ( rji_sq >= cutLJsq[itype][jtype] )
          continue;

        rji = neighbours_bonds_i[jj].r;
        rji_vec[0] = neighbours_bonds_i[jj].x;
        rji_vec[1] = neighbours_bonds_i[jj].y;
        rji_vec[2] = neighbours_bonds_i[jj].z;

        if ( rji >= cut3rebo ) {
          best = 0.0;
          testpath = false;
        } else if ( rji >= rcMax[itype][jtype] ) {
          best = 0.0;
          testpath = true;
        } else {
          best = SpRC( rji, itype, jtype );
          testpath = ( best < 1.0 );
        }

        done = false;
        if ( testpath ) {
          REBO_neighbours_num_i = REBO_neighbours_num[i];
          REBO_neighbours_list_i = REBO_neighbours_list[i];
          REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];

          for( kk = 0; ( kk < REBO_neighbours_num_i ) && !done; kk++ ) {
            k = REBO_neighbours_list_i[kk];

            if ( k == j )
              continue;

            ktype = type[k];

            rki_vec[0] = REBO_neighbours_bonds_i[kk].x;
            rki_vec[1] = REBO_neighbours_bonds_i[kk].y;
            rki_vec[2] = REBO_neighbours_bonds_i[kk].z;
            rki = REBO_neighbours_bonds_i[kk].r;
            rki_sq = REBO_neighbours_bonds_i[kk].r_sq;

            wki = ( rki_sq < rcMaxSq[itype][ktype] ) ? SpRC( rki, itype, ktype ) : 0.0;

            if ( wki <= best )
              continue;
            diff( rkj_vec, rki_vec, rji_vec );
            rkj_sq = length_sq( rkj_vec );

            if ( rkj_sq < rcMaxSq[ktype][jtype] ) {
              rkj = sqrt( rkj_sq );
              wkj = SpRC( rkj, ktype, jtype );
              if ( wki * wkj > best ) {
                best = wki * wkj;
                if ( best == 1.0 )
                  break;
              }
            }

            REBO_neighbours_num_k = REBO_neighbours_num[k];
            REBO_neighbours_list_k = REBO_neighbours_list[k];
            REBO_neighbours_bonds_k = REBO_neighbours_bonds[k];
            for( mm = 0; mm < REBO_neighbours_num_k /*&& !done*/; mm++ ) {
              m = REBO_neighbours_list_k[mm];

              if ( ( m == i ) || ( m == j ) )
                continue;

              mtype = type[m];

              rmk_vec[0] = REBO_neighbours_bonds_k[mm].x;
              rmk_vec[1] = REBO_neighbours_bonds_k[mm].y;
              rmk_vec[2] = REBO_neighbours_bonds_k[mm].z;
              rmk = REBO_neighbours_bonds_k[mm].r;
              rmk_sq = REBO_neighbours_bonds_k[mm].r_sq;

              wmk = ( rmk_sq < rcMaxSq[ktype][mtype] ) ? SpRC( rmk, ktype, mtype ) : 0.0;

              if ( wki * wmk <= best )
                continue;
              sum( rmj_vec, rmk_vec, rkj_vec );
              rmj_sq = rmj_vec[0] * rmj_vec[0] + rmj_vec[1] * rmj_vec[1] + rmj_vec[2] * rmj_vec[2];

              if ( rmj_sq < rcMaxSq[mtype][jtype] ) {
                rmj = sqrt( rmj_sq );
                wmj = SpRC( rmj, mtype, jtype );
                if ( wki * wmk * wmj > best ) {
                  best = wki * wmk * wmj;
                  if ( best == 1.0 ) {
                    done = true;
                    break;
                  }
                }
              }
            }
          }
        }

        cij = 1.0 - best;
        if ( cij == 0.0 )
          continue;

        sigwid = 0.84;
        sigcut = 3.0;
        sigmin = sigcut - sigwid;

        rljmin = sigma[itype][jtype];
        rljmax = sigcut * rljmin;
        rljmin = sigmin * rljmin;

        if ( rji > rljmax )
          slw = 0.0;
        else if ( rji > rljmin ) {
          drji = rji - rljmin;
          swidth = rljmax - rljmin;
          tee = drji / swidth;
          tee2 = tee * tee;
          slw = 1.0 - tee2 * ( 3.0 - 2.0 * tee );
        } else
          slw = 1.0;

        r2inv = 1.0 / rji_sq;
        r6inv = r2inv * r2inv * r2inv;

        vdw = r6inv * ( lj3[itype][jtype] * r6inv - lj4[itype][jtype] );
        VLJ = vdw * slw;

        Str = Sp2bLJ( rji, itype, jtype );
        VA = Str * cij * VLJ;
        if ( Str > 0.0 ) {
          scale( delscale, rji_vec, rcMin[itype][jtype] / rji );
          Stb = bond_orderLJ( i, j, delscale, rcMin[itype][jtype], rji );
        } else
          Stb = 0.0;

        eLJ += VA * Stb + ( 1.0 - Str ) * cij * VLJ;
      }
    }
    return eLJ;
  }

  double ForceField::bond_orderLJ( int i, int j, double rji_vec[3], double rji, double rji0 ) {
    int atomi, itype;
    int atomj, jtype;
    int k, atomk, ktype;
    int l, atoml, ltype;
    double wij;
    double NijC, NijH, NjiC, NjiH;
    double bij;
    double NconjtmpI, NconjtmpJ;
    double Etmp;
    double Stb;

    int REBO_neighbours_num_i;
    int *REBO_neighbours_list_i;
    vec3d *REBO_neighbours_bonds_i;
    int REBO_neighbours_num_j;
    int *REBO_neighbours_list_j;
    vec3d *REBO_neighbours_bonds_j;

    double rji_sq;
    double rki_vec[3], rki, rki_sq;
    double rlj_vec[3], rlj, rlj_sq;
    double rkj_vec[3], rkj_sq;
    double rli_vec[3], rli_sq;

    double wki, wlj;
    double cosjik, cosijl;
    double g;

    double PijS, PjiS;
    double pij, pji;

    double Nijconj;
    double piRC, Tij;

    double cos321, cos234;
    double costmp;
    double tspjik, tspijl;

    double crosskij[3];
    double crossijl[3];
    double omkijl;

    atomi = i;
    atomj = j;
    itype = type[atomi];
    jtype = type[atomj];

    wij = SpRC( rji0, itype, jtype );

    NijC = nC[atomi];
    NijH = nH[atomi];
    NjiC = nC[atomj];
    NjiH = nH[atomj];

    if ( itype == 0 )
      NjiC -= wij;
    else
      NjiH -= wij;
    if ( jtype == 0 )
      NijC -= wij;
    else
      NijH -= wij;

    bij = 0.0;
    NconjtmpI = 0.0;
    NconjtmpJ = 0.0;
    Etmp = 0.0;
    Stb = 0.0;

    REBO_neighbours_num_i = REBO_neighbours_num[i];
    REBO_neighbours_list_i = REBO_neighbours_list[i];
    REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];
    for( k = 0; k < REBO_neighbours_num_i; k++ ) {
      atomk = REBO_neighbours_list_i[k];
      if ( atomk == atomj )
        continue;
      ktype = type[atomk];

      rki_vec[0] = REBO_neighbours_bonds_i[k].x;
      rki_vec[1] = REBO_neighbours_bonds_i[k].y;
      rki_vec[2] = REBO_neighbours_bonds_i[k].z;
      rki = REBO_neighbours_bonds_i[k].r;

      wki = SpRC( rki, itype, ktype );
      cosjik = cos_theta_clamp( rji_vec, rki_vec, rji, rki );

      g = gSpline( cosjik, ( NijC + NijH ), itype );
      if ( itype == 1 )
        Etmp += wki * g * exp( 4.0 * ( rho[ktype][1] - rki - rho[jtype][1] + rji ) );
      if ( ktype == 0 )
        NconjtmpI += wki * SpN( nC[atomk] + nH[atomk] - wki );
    }

    PijS = 0.0;
    PijS = PijSpline( NijC, NijH, itype, jtype );
    pij = 1.0 / sqrt( 1.0 + Etmp + PijS );

    Etmp = 0.0;

    REBO_neighbours_num_j = REBO_neighbours_num[j];
    REBO_neighbours_list_j = REBO_neighbours_list[j];
    REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];
    for( l = 0; l < REBO_neighbours_num_j; l++ ) {
      atoml = REBO_neighbours_list_j[l];
      if ( atoml == atomi )
        continue;
      ltype = type[atoml];
      rlj_vec[0] = REBO_neighbours_bonds_j[l].x;
      rlj_vec[1] = REBO_neighbours_bonds_j[l].y;
      rlj_vec[2] = REBO_neighbours_bonds_j[l].z;
      rlj = REBO_neighbours_bonds_j[l].r;

      wlj = SpRC( rlj, jtype, ltype );
      cosijl = -cos_theta_clamp( rji_vec, rlj_vec, rji, rlj );

      g = gSpline( cosijl, NjiC + NjiH, jtype );
      if ( jtype == 1 )
        Etmp += wlj * g * exp( 4.0 * ( rho[ltype][1] - rlj - rho[itype][1] + rji ) );
      if ( ltype == 0 )
        NconjtmpJ += wlj * SpN( nC[atoml] + nH[atoml] - wlj );
    }

    PjiS = 0.0;
    PjiS = PijSpline( NjiC, NjiH, jtype, itype );
    pji = 1.0 / sqrt( 1.0 + Etmp + PjiS );

    Nijconj = 1.0 + NconjtmpI * NconjtmpI + NconjtmpJ * NconjtmpJ;
    piRC = piRCSpline( NijC + NijH, NjiC + NjiH, Nijconj, itype, jtype );
    Tij = 0.0;
    if ( itype == 0 && jtype == 0 )
      Tij = TijSpline( NijC + NijH, NjiC + NjiH, Nijconj );

    if ( fabs( Tij ) <= TOL )
      return Sp2bLJ( 0.5 * ( pij + pji ) + piRC, itype, jtype );

    Etmp = 0.0;
    REBO_neighbours_num_i = REBO_neighbours_num[i];
    REBO_neighbours_list_i = REBO_neighbours_list[i];
    REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];
    for( k = 0; k < REBO_neighbours_num_i; k++ ) {
      atomk = REBO_neighbours_list_i[k];
      ktype = type[atomk];
      if ( atomk == atomj )
        continue;
      rki_vec[0] = REBO_neighbours_bonds_i[k].x;
      rki_vec[1] = REBO_neighbours_bonds_i[k].y;
      rki_vec[2] = REBO_neighbours_bonds_i[k].z;
      rki = REBO_neighbours_bonds_i[k].r;
      cos321 = cos_theta_clamp( rji_vec, rki_vec, rji, rki );

      diff( rkj_vec, rki_vec, rji_vec );
      rkj_sq = length_sq( rkj_vec );

      rji_sq = rji * rji;
      rki_sq = rki * rki;
      costmp = 0.5 * ( rji_sq + rki_sq - rkj_sq ) / rji / rki;
      tspjik = Sp2th( costmp );

      if ( 1.0 - cos321 * cos321 <= TOL )
        continue;
      wki = SpRCP( rki, itype, ktype );

      REBO_neighbours_num_j = REBO_neighbours_num[j];
      REBO_neighbours_list_j = REBO_neighbours_list[j];
      REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];
      for( l = 0; l < REBO_neighbours_num_j; l++ ) {
        atoml = REBO_neighbours_list_j[l];
        ltype = type[atoml];
        if ( atoml == atomi || atoml == atomk )
          continue;
        rlj_vec[0] = REBO_neighbours_bonds_j[l].x;
        rlj_vec[1] = REBO_neighbours_bonds_j[l].y;
        rlj_vec[2] = REBO_neighbours_bonds_j[l].z;
        rlj = REBO_neighbours_bonds_j[l].r;
        rlj_sq = REBO_neighbours_bonds_j[l].r_sq;

        cos234 = -cos_theta_clamp( rji_vec, rlj_vec, rji, rlj );

        sum( rli_vec, rji_vec, rlj_vec );

        rli_sq = ( rli_vec[0] * rli_vec[0] ) + ( rli_vec[1] * rli_vec[1] ) + ( rli_vec[2] * rli_vec[2] );
        costmp = 0.5 * ( rji_sq + rlj_sq - rli_sq ) / rji / rlj;
        tspijl = Sp2th( costmp );

        if ( 1.0 - cos234 * cos234 <= TOL )
          continue;

        wlj = SpRC( rlj, jtype, ltype );
        cross( crosskij, rji_vec, rki_vec );
        cross( crossijl, rji_vec, rlj_vec );
        omkijl = -cos_theta_clamp( crosskij, crossijl, length( crosskij ), length( crossijl ) );
        Etmp += ( 1.0 - omkijl * omkijl ) * wki * wlj * ( 1.0 - tspjik ) * ( 1.0 - tspijl );
      }
    }

    return Sp2bLJ( 0.5 * ( pij + pji ) + piRC + Tij * Etmp, itype, jtype );
  }

}
#endif