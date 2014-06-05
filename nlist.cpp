/*
    File:   nlists.cpp
    Author: vaxquis
 */
#include <cassert>

#include "OpenREBO.h"
#include "nlist.h"

namespace OpenREBO {

  int get_uint( FILE* fp ) {
    char c;
    int i = 0;
    c = fgetc( fp );
    //assert(c!='0');
    while( isdigit( c ) ) {
      i *= 10;
      i += c - '0';
      c = fgetc( fp );
    }
    return i;
  }

  const double divs[] = {
    1E-1, 1E-2, 1E-3, 1E-4, 1E-5,
    1E-6, 1E-7, 1E-8, 1E-9, 1E-10,
    1E-11, 1E-12, 1E-13, 1E-14, 1E-15
  };

  double get_double( FILE* fp ) {
    bool sign = false;
    char c;
    double d = 0.0, d2 = 0.0;
    int div_nr = -1;
    c = fgetc( fp );
    if ( c == '-' ) {
      sign = true;
      c = fgetc( fp );
    }
    if ( c != '0' )
      while( isdigit( c ) ) {
        d *= 10;
        d += c - '0';
        c = fgetc( fp );
      }//
    else
      c = fgetc( fp );
    if ( c != '.' )
      return sign ? -d : d;
    c = fgetc( fp );
    while( isdigit( c ) ) {
      div_nr++;
      d2 *= 10;
      d2 += c - '0';
      c = fgetc( fp );
    }
    return sign ? -d - d2 * divs[div_nr] : d + d2 * divs[div_nr];
  }

#if 0
  atom** NList::fscanf_list( FILE* fp, int max_REBO_neighbors ) { // base version
    atom ** a, *at;
    double* v;
    int atom_nr, atom_type, atom_neighbor_cnt;
    fscanf( fp, "natoms %d\n", &atom_count );
    a = new atom*[atom_count];
    for( int i = 0; i < atom_count; i++ ) {
      fscanf( fp, "atom %d type %d\nnumber_of_neighbours %d\nneighbours:\n", &atom_nr, &atom_type, &atom_neighbor_cnt );
      mixed_type_flag |= atom_type;
      assert( atom_neighbor_cnt <= max_REBO_neighbours );
      at = new atom( atom_type, atom_neighbor_cnt );
      a[atom_nr] = at;
      for( int j = 0; j < atom_neighbor_cnt; j++ ) {
        bond& b = at->neighbor_bonds[j];
        //b = new bond( ); // 1000 clocks @ large
        v = b.v;
        //fscanf( fp, "%d %lf %lf %lf\n", &(b->target_nr), &v[0], &v[1], &v[2] ); // 6500 clocks ! ! !
        b.target_nr = get_uint( fp );
        for( int k = 0; k < 3; k++ )
          v[k] = get_double( fp ); // get_XXX totals @ 1500 clocks
        assert( atom_nr != b.target_nr );
        at->neighbor_count++;
      }
    }
    return a;
  }
#else
   atom** NList::fscanf_list( FILE* fp, int max_REBO_neighbors ) {
    atom ** a, *at;
    double* v;
    int atom_nr, atom_type, atom_neighbor_cnt, j;
    fscanf( fp, "natoms %d\n", &atom_count );
    a = new atom*[atom_count];
    for( int i = 0; i < atom_count; i++ ) {
      fscanf( fp, "atom %d type %d\nnumber_of_neighbours %d\nneighbours:\n", &atom_nr, &atom_type, &atom_neighbor_cnt );
      mixed_type_flag |= atom_type;
      assert( atom_neighbor_cnt <= max_REBO_neighbors );
      at = new atom( atom_type, atom_neighbor_cnt );
      a[atom_nr] = at;
      for( j = 0; j < atom_neighbor_cnt; j++ ) {
        bond& b = at->neighbor_bonds[j];
        //b = new bond( ); // 1000 clocks @ large
        v = b.v;
        //fscanf( fp, "%d %lf %lf %lf\n", &(b->target_nr), &v[0], &v[1], &v[2] ); // 6500 clocks ! ! !
        b.target_nr = get_uint( fp );
        for( int k = 0; k < 3; k++ )
          v[k] = get_double( fp ); // get_XXX totals @ 1500 clocks
        assert( atom_nr != b.target_nr );
      }
      at->neighbor_count = j;
    }
    return a;
  }
#endif
}