/*
    File:   nlists.cpp
    Author: vaxquis
 */
#include <cassert>

#include "OpenREBO.h"
#include "nlist.h"

namespace OpenREBO {

  atom** NList::read_list_from( istream& is ) {
    atom ** a, *at;
    bond* b;
    string tmp;
    int atom_nr, atom_type, atom_neighbor_count;
    is >> tmp;
    assert( tmp.compare( "natoms" ) == 0 );
    is >> atom_count;
    a = new atom*[atom_count];
    for( int i = 0; i < atom_count; i++ ) {
      is >> tmp;
      assert( tmp.compare( "atom" ) == 0 );
      is >> atom_nr;
      is >> tmp;
      assert( tmp.compare( "type" ) == 0 );
      is >> atom_type;
      is >> tmp;
      assert( tmp.compare( "number_of_neighbours" ) == 0 );
      is >> atom_neighbor_count;
      assert( atom_neighbor_count <= my_ff->max_REBO_neighbours );
      at = new atom( atom_type, atom_neighbor_count );
      a[atom_nr] = at;
      is >> tmp;
      assert( tmp.compare( "neighbours:" ) == 0 );
      for( int j = 0; j < atom_neighbor_count; j++ ) {
        b = new bond( );
        is >> b->target_nr;
        assert( atom_nr != b->target_nr );
        is >> b->v[0];
        is >> b->v[1];
        is >> b->v[2];
        b->r_sq = length_sq( b->v );
        b->r = sqrt( b->r_sq );
        b->inv_r = 1.0 / b->r;
        for( int i = 0; i < TYPE_COUNT; i++ ) {
          b->sp_rc[i] = my_ff->SpRC( b->r, atom_type, i );
          b->sp_rcP[i] = my_ff->SpRCP( b->r, atom_type, i );
        }
        at->neighbor_bonds[j] = b;
        at->neighbor_count++;
      }
    }
    return a;
  }

}