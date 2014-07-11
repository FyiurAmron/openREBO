
  struct REBO_consts {
    cl_double rcMin_CC, rcMax_CC, rcMaxP_CC, pi_div_delta_RC_CC, pi_div_delta_RCP_CC;
  }

  struct bond {
    cl_int target_nr;
    cl_double v[3]; // x, y, z;
    cl_double r, r_sq, r_inv, sp_rc[2], sp_rcP[2];
  };

  cl_double SpRC_CC( cl_double Xij, REBO_consts* s ) {
    if ( Xij >= s->rcMax_CC )
      return 0.0;
    if ( Xij <= s->rcMin_CC )
      return 1.0;
    return 0.5 + 0.5 * cos( ( Xij - s->rcMin_CC ) * s->pi_div_delta_RC_CC );
  }

  cl_double SpRCP_CC( cl_double Xij, REBO_consts* s ) {
    if ( Xij >= s->rcMaxP_CC )
      return 0.0;
    if ( Xij <= s->rcMin_CC )
      return 1.0;
    return 0.5 + 0.5 * cos( ( Xij - s->rcMin_CC ) * s->pi_div_delta_RCP_CC );
  }

__kernel void prep_neighbor( __global bond* bs, __global REBO_consts* s ) {
  bond *b = bs + get_global_id(0);

  b->r = sqrt( b->r_sq );
  b->r_inv = 1.0 / b->r;
  b->sp_rc[0] = SpRC_CC( b->r, s );
  b->sp_rcP[0] = SpRCP_CC( b->r, s );
}