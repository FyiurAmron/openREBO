/*
    File:   nlists.h
    Author: vaxquis
 */
#pragma once

#ifndef NLIST_H
#define	NLIST_H

#include <istream>
#include <fstream>
#include <string>

using std::string;
using std::ifstream;

#include "const.h"
#include "vector.h"

namespace OpenREBO {
  using std::istream;

  class AIREBO;

  struct bond {
    int target_nr;
    double v[3]; // x, y, z;
    double r, r_sq, r_inv, sp_rc[TYPE_COUNT], sp_rcP[TYPE_COUNT];
  };

  class neighbor_tracker {
  public:
    int neighbor_count;
    bond** neighbor_bonds;

    neighbor_tracker( int max_neighbor_count ) {
      neighbor_bonds = new bond*[max_neighbor_count];
      this->neighbor_count = 0;
    }

    virtual ~neighbor_tracker( ) {
      delete [] neighbor_bonds;
    }

  };

  class atom {
  public:
    int neighbor_count;
    bond* neighbor_bonds;

    int type;
    double nC, nH, nTotal;

    atom( int atom_type, int max_neighbor_count ) {
      neighbor_bonds = new bond[max_neighbor_count];
      this->neighbor_count = 0;
      assert( atom_type == 0 || atom_type == 1 ); // note: 0 for C, 1 for H *by contract*; -1 (LAMMPS NULL) disallowed
      type = atom_type;
      nC = 0.0;
      nH = 0.0;
      nTotal = 0.0;
    }

    virtual ~atom( ) {
      /*
      for( int i = 0; i < neighbor_count; i++ )
        delete neighbor_bonds[i];
       */
      delete [] neighbor_bonds;
    }

    atom& operator =(const atom& right ) = delete;
  };

  class NList {
  public:
    atom** atoms;
    int atom_count;
    int mixed_type_flag;

    NList( const string& filename, int max_REBO_neighbors ) {
      mixed_type_flag = 0;
      //ifstream ifs( filename );
      //atoms = read_list_from( ifs );
      FILE* fp = fopen( filename.c_str( ), "r" );
      atoms = fscanf_list( fp, max_REBO_neighbors );
      fclose( fp );
    }

    virtual ~NList( ) {
      for( int i = 0; i < atom_count; i++ )
        delete atoms[i];
      delete [] atoms;
    }

    bool is_pure_C( ) {
      return mixed_type_flag == 0;
    }

  private:
    atom** read_list_from_stream( istream& is );
    atom** fscanf_list( FILE* fp, int max_REBO_neighbors );
  };

  class RNList {
  public:
    neighbor_tracker** neighbors;
    int atom_count;

    RNList( int atom_count, int max_neighbors ) {
      this->atom_count = atom_count;
      neighbors = new neighbor_tracker*[atom_count];
      for( int i = 0; i < atom_count; i++ )
        neighbors[i] = new neighbor_tracker( max_neighbors );
    }

    virtual ~RNList( ) {
      for( int i = 0; i < atom_count; i++ )
        delete neighbors[i]; // don't do it since it only stores *existing* bonds created for 'atom' by NList loader
      delete [] neighbors;
    }
  };

}

#endif	/* NLIST_H */

