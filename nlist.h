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

namespace AIREBO {
  using std::istream;

  class ForceField;

  /*
  typedef struct {
    double x, y, z;
    double r, r_sq;
  } vec3d;
   */
  class bond {
  public:
    int target_nr;
    double v[3]; // x, y, z;
    double r, r_sq, sp_rc[TYPE_COUNT];
  };

  class neighbor_tracker {
  public:
    int neighbor_count;//, *neighbor_list;
    bond** neighbor_bonds;

    virtual ~neighbor_tracker( ) {
      //delete [] neighbor_list;
      for( int i = 0; i < neighbor_count; i++ )
        delete neighbor_bonds[i];
      delete [] neighbor_bonds;
    }

  };

  class atom : public neighbor_tracker {
  public:
    int type; // note: 0 for C, 1 for H *by contract*; -1 (LAMMPS NULL) disallowed
    double nC, nH, nTotal;

    atom( ) {
      nC = 0.0;
      nH = 0.0;
      nTotal = 0.0;
    }

    atom& operator =(const atom& right ) = delete;
  };

  class NList {
  public:
    atom** atoms;
    int atom_count;
    const ForceField* my_ff;

    NList( const string& filename, const ForceField* my_ff ) {
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
      neighbor_tracker *nt;
      for( int i = 0; i < atom_count; i++ ) {
        nt = new neighbor_tracker( );
        nt->neighbor_count = 0;
        //nt->neighbor_list = new int[max_neighbors];
        nt->neighbor_bonds = new bond*[max_neighbors];
        neighbors[i] = nt;
      }
    }

    virtual ~RNList( ) {
      for( int i = 0; i < atom_count; i++ )
        delete neighbors[i];
    }
  };

}

#endif	/* NLIST_H */

