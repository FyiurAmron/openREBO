// based on LAMMPS implementation of AIREBO force field

#include "airebo_force_field.h"

namespace AIREBO {

  void ForceField::REBO_neighbours( ) {
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

  void ForceField::E_REBO( ) {
    int i_type, j_type;
    double r_ij, w_ij, VA;
    double *Q_i, *A_i, *alpha_i;

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

        energy_rebo += w_ij * ( 1.0 + Q_i[j_type] / r_ij ) * A_i[j_type] * exp( -alpha_i[j_type] * r_ij )
                + bond_order( i, j, b ) * VA;
      }
    }
  }

  double ForceField::bond_order( int a1, int a2, const bond* b ) { //const double r_ij_vec[3], double r_ij ) {
    int type1, type2;
    int i_type;
    int a_l, l_type;

    /*
    int REBO_neighbours_num_i;
    int *REBO_neighbours_list_i;
    vec3d *REBO_neighbours_bonds_i;
    int REBO_neighbours_num_j;
    int *REBO_neighbours_list_j;
    vec3d *REBO_neighbours_bonds_j;
     */
    double *r1i_vec, *r12_vec, r12, r1i;
    /*
        double rij_vec[3], rij, rij_sq;
        double rki_vec[3], rki, rki_sq;
        double rlj_vec[3], rlj, rlj_sq;
        double rkj_vec[3], rkj_sq;
        double rli_vec[3], rli_sq;
     */
    double w12, w1i, wlj;

    double g;

    double NijC, NijH, NjiC, NjiH, NijTotal, NjiTotal;

    double NconjtmpI, NconjtmpJ, Nijconj;

    double pij, pji;
    double piRC, Tij, Etmp;

    double cos321, sin321;
    double cos234, sin234;

    double tspjik, tspijl;
    double cross321[3], cross234[3];

    double cwnum, cwnom;
    double om1234;

    atom *atom_i = atom_list->atoms[atom_i],
            *atom_j = atom_list->atoms[atom_j];
    type1 = atom_i->type;
    type2 = atom_j->type;

    r12_vec = b->v;
    r12 = b->r;
    w12 = b->sp_rc[type2];

    NijC = atom_i->nC;
    NijH = atom_i->nH;
    NjiC = atom_j->nC;
    NjiH = atom_j->nH;

    if ( type1 == 0 )
      NjiC -= w12;
    else
      NjiH -= w12;

    if ( type2 == 0 )
      NijC -= w12;
    else
      NijH -= w12;

    NijTotal = NijC + NijH;
    NjiTotal = NjiC + NjiH;

    NconjtmpI = 0.0;
    NconjtmpJ = 0.0;
    Etmp = 0.0;

    neighbor_tracker *nt = rebo_atom_list->neighbors[a1];
    bond* b2;
    for( int i = 0, max_i = nt->neighbor_count; i < max_i; i++ ) {
      b2 = nt->neighbor_bonds[i];
      atom_i = b2->target_nr;
      if ( atom_i == a2 )
        continue;

      i_type = atom_list->atoms[atom_i]->type;
      r1i_vec = b2->v;
      r1i = b2->r;
      w1i = b2->sp_rc[i_type];

      Etmp += gSpline( cos_theta_clamp( r12_vec, r1i_vec, r12, r1i ), NijTotal, type1 )
              * ( ( type1 == 1 ) ?
              ( w1i * exp( 4.0 * ( ( i_type != type2 ) ? ( rho[i_type][1] - rho[type2][1] + r12 - r1i ) : r12 - r1i ) ) )
              : w1i );

      if ( i_type == 0 )
        NconjtmpI += w1i * SpN( atom_i->nTotal - w1i );
    }
    ////////////////////////////////////////////////////////
    pij = pow( 1.0 + Etmp + PijSpline( NijC, NijH, type1, type2 ), -0.5 );

    Etmp = 0.0;

    REBO_neighbours_num_j = REBO_neighbours_num[j];
    REBO_neighbours_list_j = REBO_neighbours_list[j];
    REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];

    for( int j = 0; j < REBO_neighbours_num_j; j++ ) {
      a_l = REBO_neighbours_list_j[j];
      if ( a_l == a1 )
        continue;

      l_type = type[a_l];
      rlj_vec[0] = REBO_neighbours_bonds_j[j].x;
      rlj_vec[1] = REBO_neighbours_bonds_j[j].y;
      rlj_vec[2] = REBO_neighbours_bonds_j[j].z;
      rlj = REBO_neighbours_bonds_j[j].r;

      wlj = SpRC( rlj, type2, l_type );

      g = gSpline( cos_theta_clamp( r_ij_vec, rlj_vec, r12, rlj ), NjiTotal, type2 );
      if ( type2 == 1 )
        Etmp += wlj * g * exp( 4.0 * ( ( rho[l_type][1] - rlj ) - ( rho[type1][1] - r12 ) ) );
      ERROR see above
      if ( l_type == 0 )
        NconjtmpJ += wlj * SpN( nC[a_l] + nH[a_l] - wlj );
    }

    pji = pow( 1.0 + Etmp + PijSpline( NjiC, NjiH, type2, type1 ), -0.5 );

    Nijconj = 1.0 + NconjtmpI * NconjtmpI + NconjtmpJ * NconjtmpJ;
    piRC = piRCSpline( NijTotal, NjiTotal, Nijconj, type1, type2 );

    Tij = 0.0;
    if ( ( type1 == 0 ) && ( type2 == 0 ) )
      Tij = TijSpline( NijTotal, NjiTotal, Nijconj );
    Etmp = 0.0;

    if ( fabs( Tij ) <= TOL )
      return 0.5 * ( pij + pji ) + piRC;

    neg( rij_vec, r_ij_vec );
    rij = r12;
    rij_sq = rij * rij;

    REBO_neighbours_num_i = REBO_neighbours_num[i];
    REBO_neighbours_list_i = REBO_neighbours_list[i];
    REBO_neighbours_bonds_i = REBO_neighbours_bonds[i];
    for( int i = 0; i < REBO_neighbours_num_i; i++ ) {
      atom_i = REBO_neighbours_list_i[i];
      i_type = type[atom_i];
      if ( atom_i == a2 )
        continue;
      rki_vec[0] = REBO_neighbours_bonds_i[i].x;
      rki_vec[1] = REBO_neighbours_bonds_i[i].y;
      rki_vec[2] = REBO_neighbours_bonds_i[i].z;
      rki = REBO_neighbours_bonds_i[i].r;
      rki_sq = REBO_neighbours_bonds_i[i].r_sq;

      cos321 = cos_theta_clamp( rki_vec, rij_vec, rki, rij );
      sin321 = cos321 * cos321;
      if ( sin321 == 1.0 )
        continue;
      sin321 = sqrt( 1.0 - sin321 );

      diff( rkj_vec, rki_vec, r_ij_vec );
      rkj_sq = length_sq( rkj_vec );
      w1i = SpRCP( rki, type1, i_type );
      tspjik = Sp2th( 0.5 * ( rij_sq + rki_sq - rkj_sq ) / ( rij * rki ) );

      REBO_neighbours_num_j = REBO_neighbours_num[j];
      REBO_neighbours_list_j = REBO_neighbours_list[j];
      REBO_neighbours_bonds_j = REBO_neighbours_bonds[j];
      for( int j = 0; j < REBO_neighbours_num_j; j++ ) {
        a_l = REBO_neighbours_list_j[j];
        l_type = type[a_l];
        if ( ( a_l == a1 ) || ( a_l == atom_i ) )
          continue;
        rlj_vec[0] = REBO_neighbours_bonds_j[j].x;
        rlj_vec[1] = REBO_neighbours_bonds_j[j].y;
        rlj_vec[2] = REBO_neighbours_bonds_j[j].z;
        rlj = REBO_neighbours_bonds_j[j].r;
        rlj_sq = REBO_neighbours_bonds_j[j].r_sq;

        cos234 = cos_theta_clamp( rij_vec, rlj_vec, rij, rlj );
        sin234 = cos234 * cos234;
        if ( sin234 == 1.0 )
          continue;
        sin234 = sqrt( 1.0 - sin234 );

        wlj = SpRCP( rlj, type2, l_type );
        sum( rli_vec, r_ij_vec, rlj_vec );
        rli_sq = length_sq( rli_vec );

        tspijl = Sp2th( 0.5 * ( rij_sq + rlj_sq - rli_sq ) / ( rij * rlj ) );

        cross( cross321, rij_vec, rki_vec );
        cross( cross234, r_ij_vec, rlj_vec );
        cwnum = dot( cross321, cross234 );
        cwnom = rki * rlj * r12 * r12 * sin321 * sin234;
        om1234 = cwnum / cwnom;

        Etmp += ( ( 1.0 - om1234 * om1234 ) * w1i * wlj ) * ( 1.0 - tspjik ) * ( 1.0 - tspijl );
      }
    }

    return 0.5 * ( pij + pji ) + piRC + Tij * Etmp;
  } // bond_order
}
