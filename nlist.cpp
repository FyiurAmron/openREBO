/*
    File:   nlists.cpp
    Author: vaxquis
 */
#include <cassert>

#include "airebo_force_field.h"
#include "nlist.h"

namespace AIREBO {

  atom** NList::read_list_from( istream& is ) {
    atom ** a, *at;
    bond* b;
    string tmp;
    int atom_nr;
    is >> tmp;
    assert( tmp.compare( "natoms" ) == 0 );
    is >> atom_count;
    a = new atom*[atom_count];
    for( int i = 0; i < atom_count; i++ ) {
      at = new atom( );
      is >> tmp;
      assert( tmp.compare( "atom" ) == 0 );
      is >> atom_nr;
      is >> tmp;
      assert( tmp.compare( "type" ) == 0 );
      is >> at->type;
      is >> tmp;
      assert( tmp.compare( "number_of_neighbours" ) == 0 );
      is >> at->neighbor_count;
      assert( at->neighbor_count <= my_ff->max_REBO_neighbours );
      a[atom_nr] = at;
      at->neighbor_bonds = new bond*[at->neighbor_count];
      is >> tmp;
      assert( tmp.compare( "neighbours:" ) == 0 );
      for( int j = 0; j < at->neighbor_count; j++ ) {
        b = new bond( );
        is >> b->target_nr;
        assert( atom_nr != b->target_nr );
        is >> b->v[0];
        is >> b->v[1];
        is >> b->v[2];
        b->r_sq = length( b->v );
        b->r = sqrt( b->r_sq );
        for( int i = 0; i < TYPE_COUNT; i++ )
          b->sp_rc[i] = my_ff->SpRC( b->r, at->type, i );
        at->neighbor_bonds[j] = b;
      }
    }
    return a;
  }

}