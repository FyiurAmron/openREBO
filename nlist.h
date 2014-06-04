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

  class bond {
  public:
    int target_nr;
    double v[3]; // x, y, z;
    double r, r_sq, inv_r, sp_rc[TYPE_COUNT], sp_rcP[TYPE_COUNT];
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

  class atom : public neighbor_tracker {
  public:
    int type;
    double nC, nH, nTotal;

    atom( int atom_type, int neighbor_count ) : neighbor_tracker( neighbor_count ) {
      assert( atom_type == 0 || atom_type == 1 ); // note: 0 for C, 1 for H *by contract*; -1 (LAMMPS NULL) disallowed
      type = atom_type;
      nC = 0.0;
      nH = 0.0;
      nTotal = 0.0;
    }

    virtual ~atom( ) {
      for( int i = 0; i < neighbor_count; i++ )
        delete neighbor_bonds[i];
    }

    atom& operator =(const atom& right ) = delete;
  };

  class NList {
  public:
    atom** atoms;
    int atom_count;
    const AIREBO* my_ff;

    NList( const string& filename, const AIREBO* my_ff ) {
      this->my_ff = my_ff;
      ifstream ifs( filename );
      atoms = read_list_from( ifs );
    }

    /*
    NList( istream& is ) {
      atoms = read_list_from( is );
    }
     */
    virtual ~NList( ) {
      for( int i = 0; i < atom_count; i++ )
        delete atoms[i];
      delete [] atoms;
    }

  private:
    atom** read_list_from( istream& is );
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

