  struct bond {
    cl_int target_nr;
    cl_double v[3]; // x, y, z;
    cl_double r, r_sq, r_inv, sp_rc[TYPE_COUNT], sp_rcP[TYPE_COUNT];
  };

cl_double length_sq( cl_double* v ) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

cl_double Sp( cl_double Xij, cl_double Xmin, cl_double Xmax )  {
    if ( t <= Xmin )
        return 1.0;
    else if ( t >= Xmax )
        return -1.0;
    return 0.5 * ( 1.0 + cos( ( Xij - Xmin ) / ( Xmax - Xmin ) * PI ) );
}

__kernel void prep_neighbor( __global bond* bs, const cl_double rcMaxSq,
  const cl_double rcMin, const cl_double rcMax, const cl_double rcPMax  ) { 
  bond *b = bs + get_global_id(0);
  cl_double r = length_sq( b->v );
  if ( r >= rcMaxSq )
    return;

  //nt_i_bonds[tmp] = &b;
  b->r_sq = r;
  b->r = half_sqrt( b->r_sq );
  b->r_inv = 1.0 / b->r;
  //b->sp_rc[0] = Sp( b->r, rcMin, rcMax ); // etc
  //b->sp_rc[0] = Sp( b->r, rcMin, rcMaxP ); // etc

  atom_i->nC += b.sp_rc[0];
  tmp++;
}