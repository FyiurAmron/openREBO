/*
   File:   E_REBO_C.cpp
   Author: vaxquis
   based on LAMMPS implementation of AIREBO potential
   Carbon-only optimized calculations
 */

#include "OpenREBO.h"

namespace OpenREBO {

  void AIREBO::REBO_neighbours_C( ) {
    rebo_atom_list = new RNList( atom_list->atom_count, max_REBO_neighbours );

    int tmp;
    double r;

    atom *atom_i;
    neighbor_tracker* nt_i;
    bond **nt_i_bonds;

    double rcMaxSqIJ = rcMaxSq[0][0];
    for( int i = 0, max = atom_list->atom_count; i < max; i++ ) {
      atom_i = atom_list->atoms[i];
      nt_i = rebo_atom_list->neighbors[i];
      nt_i_bonds = nt_i->neighbor_bonds;
      tmp = 0;
      for( int j = 0, max = atom_i->neighbor_count; j < max; j++ ) {
        bond &b = atom_i->neighbor_bonds[j];
        r = length_sq( b.v );
        if ( r >= rcMaxSqIJ )
          continue;

        nt_i_bonds[tmp] = &b;
        b.r_sq = r;
        b.r = sqrt( b.r_sq );
        b.r_inv = 1.0 / b.r;
        b.sp_rc[0] = SpRC_CC( b.r );
        b.sp_rcP[0] = SpRCP_CC( b.r );

        atom_i->nC += b.sp_rc[0];
        tmp++;
      }
      atom_i->nTotal = atom_i->nC;
      nt_i->neighbor_count = tmp;
    }
  }

  double AIREBO::E_REBO_C( ) {
    int j_nr;
    double r_ij, w_ij, VA;
    double eREBO = 0.0;

    neighbor_tracker *nt_i;
    bond **nt_i_bonds, *b;

    double Q00 = Q[0][0],
            A00 = A[0][0],
            alpha00 = -alpha[0][0];
    for( int i = 0, max = atom_list->atom_count; i < max; i++ ) {
      nt_i = rebo_atom_list->neighbors[i];
      nt_i_bonds = nt_i->neighbor_bonds;
      for( int j = 0, max = nt_i->neighbor_count; j < max; j++ ) {
        b = nt_i_bonds[j];
        j_nr = b->target_nr;
        if ( i > j_nr )
          continue;

        w_ij = nt_i_bonds[j]->sp_rc[0];
        if ( w_ij <= TOL )
          continue;

        r_ij = b->r;
        VA = 0.0;
        for( int k = 0; k < 3; k++ )
          VA -= w_ij * BIJc[0][0][k] * exp( -beta[0][0][k] * r_ij );

        eREBO += w_ij * ( 1.0 + Q00 * b->r_inv ) * A00 * exp( alpha00 * r_ij )
                + bond_order_C( i, j_nr, b ) * VA;
      }
    }
    return eREBO;
  }

  double AIREBO::bond_order_1_C( int a1, int a2, const double r12_vec[3], double r12, double r12_inv,
          double NC, double &p, const neighbor_tracker * nt ) {
    bond* b2;
    atom *atom_i;

    int a_i;
    double Nconj_tmp = 0.0, w_bi;
    double Etmp = 1.0 + PijSplineC_only( NC );

    for( int i = 0, max_i = nt->neighbor_count; i < max_i; i++ ) {
      b2 = nt->neighbor_bonds[i];
      a_i = b2->target_nr;
      if ( a_i == a2 )
        continue;

      atom_i = atom_list->atoms[a_i];
      w_bi = b2->sp_rc[0];
      Etmp += gSplineC( cos_theta( r12_vec, b2->v, r12_inv * b2->r_inv ), NC ) * w_bi;

      Nconj_tmp += w_bi * SpN( atom_i->nTotal - w_bi );
    }
    p = 1.0 / sqrt( Etmp );
    return Nconj_tmp;
  }

  double AIREBO::bond_order_C( int a1, int a2, const bond * b ) {
    atom *atom1 = atom_list->atoms[a1],
            *atom2 = atom_list->atoms[a2];

    const double *r12_vec = b->v;
    double r21_vec[] = { -r12_vec[0], -r12_vec[1], -r12_vec[2] };
    double r12 = b->r,
            r12_sq = b->r_sq,
            r12_inv = b->r_inv;

    double N12C = atom1->nC - b->sp_rc[0],
            N21C = atom2->nC - b->sp_rc[0];

    neighbor_tracker
            *nt1 = rebo_atom_list->neighbors[a1],
            *nt2 = rebo_atom_list->neighbors[a2];

    double p12 = 0.0, p21 = 0.0;
    double Nconj_tmp1 = bond_order_1_C( a1, a2, r12_vec, r12, r12_inv, N12C, p12, nt1 ),
            Nconj_tmp2 = bond_order_1_C( a2, a1, r21_vec, r12, r12_inv, N21C, p21, nt2 );
    double N12conj = 1.0 + Nconj_tmp1 * Nconj_tmp1 + Nconj_tmp2 * Nconj_tmp2;
    double ret = piRCSplineCC( N12C, N21C, N12conj ) + 0.5 * ( p12 + p21 );

    double T12 = TijSpline( N12C, N21C, N12conj );
    if ( fabs( T12 ) <= TOL )
      return ret;

    bond *b2, *b3;
    double Etmp = 0.0, EwTmp,
            b2_r_inv, b3_r_inv, r12b2_inv, r12b3_inv,
            dot321, dot234, cos321, sin321, cos234, sin234, om1234;
    double rTmp_vec[3], *b2_v;
    int a_i, a_j;
    int max_j = nt2->neighbor_count;
    for( int i = 0, max_i = nt1->neighbor_count; i < max_i; i++ ) {
      b2 = nt1->neighbor_bonds[i];
      a_i = b2->target_nr;
      if ( a_i == a2 )
        continue;
      b2_r_inv = b2->r_inv;
      r12b2_inv = b2_r_inv * r12_inv;
      b2_v = b2->v;

      dot321 = dot( b2_v, r12_vec );
      cos321 = dot321 * r12b2_inv;
      sin321 = cos321 * cos321;
      if ( sin321 == 1.0 )
        continue;
      sin321 = sqrt( 1.0 - sin321 );

      sum( rTmp_vec, b2_v, r21_vec );
      EwTmp = b2->sp_rcP[0]
              * ( 1.0 - Sp2th( 0.5 * ( r12_sq + b2->r_sq - length_sq( rTmp_vec ) ) * r12b2_inv ) );
      for( int j = 0; j < max_j; j++ ) {
        b3 = nt2->neighbor_bonds[j];
        b3_r_inv = b3->r_inv;
        r12b3_inv = r12_inv * b3_r_inv;

        a_j = b3->target_nr;

        if ( a_j == a1 || a_j == a_i )
          continue;

        dot234 = dot( b3->v, r21_vec );
        cos234 = dot234 * r12b3_inv;
        sin234 = cos234 * cos234;
        if ( sin234 == 1.0 )
          continue;
        sin234 = sqrt( 1.0 - sin234 );

        double cross321[3], cross234[3];
        cross( cross321, r21_vec, b2_v );
        cross( cross234, b3->v, r21_vec );
        om1234 = ( dot234 * dot321 * r12_inv * r12_inv + dot( b2_v, b3->v ) ) * b2_r_inv * b3_r_inv / ( sin321 * sin234 );

        sum( rTmp_vec, r12_vec, b3->v );
        Etmp += EwTmp * ( 1.0 - om1234 * om1234 ) * b3->sp_rcP[0]
                * ( 1.0 - Sp2th( 0.5 * ( r12_sq + b3->r_sq - length_sq( rTmp_vec ) ) * r12b3_inv ) );
      }
    }

    return ret + T12 * Etmp;
  } // bond_order
}
