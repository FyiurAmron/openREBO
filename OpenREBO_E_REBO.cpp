/*
    File:   OpenREBO_E_REBO.cpp
    Author: vaxquis
    based on LAMMPS implementation of AIREBO potential
 */

#include "OpenREBO.h"

namespace OpenREBO {

  void AIREBO::REBO_neighbours( ) {
    rebo_atom_list = new RNList( atom_list->atom_count, max_REBO_neighbours );
    double* rcMaxSqI;
    int j_type, tmp;

    atom *atom_i;
    neighbor_tracker* nt_i;
    bond** nt_i_bonds;

    for( int i = 0, max = atom_list->atom_count; i < max; i++ ) {
      atom_i = atom_list->atoms[i];
      nt_i = rebo_atom_list->neighbors[i];
      nt_i_bonds = nt_i->neighbor_bonds;
      rcMaxSqI = rcMaxSq[atom_i->type];
      tmp = 0;
      for( int j = 0, max = atom_i->neighbor_count; j < max; j++ ) {
        bond* bond_ij = atom_i->neighbor_bonds[j];
        j_type = atom_list->atoms[bond_ij->target_nr]->type;

        if ( bond_ij->r_sq >= rcMaxSqI[j_type] )
          continue;

        nt_i_bonds[tmp] = bond_ij;

        if ( j_type == 0 ) // C
          atom_i->nC += bond_ij->sp_rc[0];
        else // H
          atom_i->nH += bond_ij->sp_rc[1];
        tmp++;
      }
      atom_i->nTotal = atom_i->nC + atom_i->nH;
      nt_i->neighbor_count = tmp;
    }
  }

  double AIREBO::E_REBO( ) {
    int i_type, j_type;
    double r_ij, w_ij, VA;
    double *Q_i, *A_i, *alpha_i;
    double eREBO = 0.0;

    atom *atom_j;
    neighbor_tracker *nt_i;
    bond **nt_i_bonds, *b;

    for( int i = 0, max = atom_list->atom_count; i < max; i++ ) {
      i_type = atom_list->atoms[i]->type;
      nt_i = rebo_atom_list->neighbors[i];
      nt_i_bonds = nt_i->neighbor_bonds;
      Q_i = Q[i_type];
      A_i = A[i_type];
      alpha_i = alpha[i_type];
      //beta_i = beta[i_type];
      //BIJc_i = BIJc[i_type];
      for( int j = 0, max = nt_i->neighbor_count; j < max; j++ ) {
        b = nt_i_bonds[j];
        atom_j = atom_list->atoms[b->target_nr];
        j_type = atom_j->type;

        w_ij = nt_i_bonds[j]->sp_rc[j_type];
        if ( w_ij <= TOL )
          continue;

        r_ij = b->r;
        VA = 0.0;
        for( int k = 0; k < 3; k++ )
          VA -= w_ij * BIJc[i_type][j_type][k] * exp( -beta[i_type][j_type][k] * r_ij );

        eREBO += w_ij * ( 1.0 + Q_i[j_type] * b->inv_r ) * A_i[j_type] * exp( -alpha_i[j_type] * r_ij )
                + bond_order( i, j, b ) * VA;
      }
    }
    return eREBO;
  }

  double AIREBO::bond_order_1( int a1, int a2, int type1, int type2, const double r12_vec[3], double r12, double inv_r12,
          double NC, double NH, double NTotal, double &p, const neighbor_tracker *nt ) {
    bond* b2;
    atom *atom_i;

    int a_i, i_type;
    double Nconj_tmp = 0.0, rho_bi, w_bi;
    double Etmp = 1.0 + PijSpline( NC, NH, type1, type2 );

    for( int i = 0, max_i = nt->neighbor_count; i < max_i; i++ ) {
      b2 = nt->neighbor_bonds[i];
      a_i = b2->target_nr;
      if ( a_i == a2 )
        continue;

      atom_i = atom_list->atoms[a_i];
      i_type = atom_i->type;
      rho_bi = r12 - b2->r;
      w_bi = b2->sp_rc[i_type];

      Etmp += gSpline( cos_theta_clamp( r12_vec, b2->v, inv_r12 * b2->inv_r ), NTotal, type1 )
              * ( ( type1 == 1 ) ?
              ( w_bi * exp( 4.0 * ( ( i_type != type2 ) ? ( rho[i_type][1] - rho[type2][1] + rho_bi ) : rho_bi ) ) )
              : w_bi );

      if ( i_type == 0 )
        Nconj_tmp += w_bi * SpN( atom_i->nTotal - w_bi );
    }
    p = 1.0 / sqrt( Etmp );
    return Nconj_tmp;
  }

  double AIREBO::bond_order( int a1, int a2, const bond* b ) {
    atom *atom1 = atom_list->atoms[a1],
            *atom2 = atom_list->atoms[a2];
    int type1 = atom1->type,
            type2 = atom2->type;

    const double *r12_vec = b->v;
    double r21_vec[] = { -r12_vec[0], -r12_vec[1], -r12_vec[2] };
    double r12 = b->r,
            r12_sq = b->r_sq,
            inv_r12 = b->inv_r;

    double N12C = atom1->nC,
            N12H = atom1->nH,
            N21C = atom2->nC,
            N21H = atom2->nH;

    double w12 = b->sp_rc[type2];
    if ( type1 == 0 )
      N21C -= w12;
    else
      N21H -= w12;

    if ( type2 == 0 )
      N12C -= w12;
    else
      N12H -= w12;

    double N12Total = N12C + N12H,
            N21Total = N21C + N21H;

    neighbor_tracker *nt1 = rebo_atom_list->neighbors[a1];
    neighbor_tracker *nt2 = rebo_atom_list->neighbors[a2];

    double p12 = 0.0, p21 = 0.0;
    double Nconj_tmp1 = bond_order_1( a1, a2, type1, type2, r12_vec, r12, inv_r12, N12C, N12H, N12Total, p12, nt1 ),
            Nconj_tmp2 = bond_order_1( a2, a1, type2, type1, r21_vec, r12, inv_r12, N21C, N21H, N21Total, p21, nt2 );
    double N12conj = 1.0 + Nconj_tmp1 * Nconj_tmp1 + Nconj_tmp2 * Nconj_tmp2;
    double ret = piRCSpline( N12Total, N21Total, N12conj, type1, type2 ) + 0.5 * ( p12 + p21 );

    if ( type1 != 0 || type2 != 0 )
      return ret;

    double T12 = TijSpline( N12Total, N21Total, N12conj );
    if ( fabs( T12 ) <= TOL )
      return ret;

    bond *b2, *b3;
    //atom *atom_i;
    int a_i, a_j;
    double Etmp = 0.0, EwTmp, b2_inv_r;
    double cos321, sin321, cos234, sin234, om1234;
    double rTmp_vec[3], cross321[3], cross234[3], *b2_v;
    int max_j = nt2->neighbor_count;
    for( int i = 0, max_i = nt1->neighbor_count; i < max_i; i++ ) {
      b2 = nt1->neighbor_bonds[i];
      a_i = b2->target_nr;
      if ( a_i == a2 )
        continue;
      b2_inv_r = b2->inv_r;
      b2_v = b2->v;

      cos321 = cos_theta_clamp( b2_v, r12_vec, b2_inv_r * inv_r12 );
      sin321 = cos321 * cos321;
      if ( sin321 == 1.0 )
        continue;
      sin321 = sqrt( 1.0 - sin321 );

      sum( rTmp_vec, b2_v, r21_vec );
      EwTmp = b2->sp_rcP[atom_list->atoms[a_i]->type]
              * ( 1.0 - Sp2th( 0.5 * ( r12_sq + b2->r_sq - length_sq( rTmp_vec ) ) * inv_r12 * b2->inv_r ) );
      for( int j = 0; j < max_j; j++ ) {
        b3 = nt2->neighbor_bonds[j];
        a_j = b3->target_nr;

        if ( a_j == a1 || a_j == a_i )
          continue;

        cos234 = cos_theta_clamp( b3->v, r21_vec, inv_r12 * b3->inv_r );
        sin234 = cos234 * cos234;
        if ( sin234 == 1.0 )
          continue;
        sin234 = sqrt( 1.0 - sin234 );


        cross( cross321, r21_vec, b2_v ); // simplify by linear alg. properties?
        cross( cross234, r12_vec, b3->v ); // ditto
        om1234 = dot( cross321, cross234 ) * b2_inv_r * b3->inv_r / ( r12_sq * sin321 * sin234 ); // ditto

        sum( rTmp_vec, r12_vec, b3->v );
        Etmp += EwTmp * ( 1.0 - om1234 * om1234 ) * b3->sp_rcP[atom_list->atoms[a_j]->type]
                * ( 1.0 - Sp2th( 0.5 * ( r12_sq + b3->r_sq - length_sq( rTmp_vec ) ) * inv_r12 * b3->inv_r ) );
      }
    }

    return ret + T12 * Etmp;
  } // bond_order
}
