#pragma once

#ifndef AIREBO_FORCE_FIELD_H
#define AIREBO_FORCE_FIELD_H

/*
 *   written by Szymon Winczewski
 *   based on LAMMPS implementation of AIREBO force field
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#define AIREBO_DEBUG

const double PI = M_PI,
        TOL = 1.0e-9;
const int TYPE_COUNT = 2;

typedef struct {
  double x, y, z;
  double r, r_sq;
} vec3d;

using std::ifstream;
using std::getline;
using std::string;
using std::endl;
using std::istringstream;
using std::cout;

class AIREBOForceField {
public:
  AIREBOForceField( string file_name, double cutlj,
                    bool ljflag, bool torflag, int max_number_of_REBO_neighbours );
  ~AIREBOForceField( );

  double getCutoffRadius( );
  double compute( int number_of_atoms, int *type,
                  int *neighbours_num, int **neighbours_list,
                  vec3d **neighbours_bonds );

  double getTotalEnergy( );
  double getREBOEnergy( );
  double getLJEnergy( );
  double getTORSIONEnergy( );

private:
  // ogolne parametry obliczen
  int natoms;
  int *type;

  // sumy funkcji wagowych (czlon REBO)
  double *nC, *nH;

  // sasiedzi REBO
  int *REBO_neighbours_num;
  int **REBO_neighbours_list;
  vec3d **REBO_neighbours_bonds;

  // wszyscy sasiedzi
  int *neighbours_num;
  int **neighbours_list;
  vec3d **neighbours_bonds;

  // promien odciecia dla czlonu LJ, wyrazony w wielokrotnosciach parametru sigma potencjalu LJ
  double cutlj;
  double **cutljsq;
  double **lj1, **lj2, **lj3, **lj4;
  double cut3rebo, cutljrebo, cutljrebosq;
  double cutmax;

  // flagi: czy uwzgledniac czlon lj oraz czlon torsyjny
  bool ljflag, torflag;
  int max_number_of_REBO_neighbours;

  double total_energy, energy_rebo, energy_lj, energy_torsion;

  // parametry potencjalu AIREBO
  // a) czlon REBO
  double rcMin[TYPE_COUNT][TYPE_COUNT], rcMax[TYPE_COUNT][TYPE_COUNT],
  rcMaxSq[TYPE_COUNT][TYPE_COUNT], rcMaxP[TYPE_COUNT][TYPE_COUNT],
  pi_div_delta_rc[TYPE_COUNT][TYPE_COUNT], pi_div_delta_rcp[TYPE_COUNT][TYPE_COUNT],
  pi_div_delta_N[TYPE_COUNT][TYPE_COUNT], pi_div_delta_NC[TYPE_COUNT][TYPE_COUNT];

  double smin;
  double Nmin, Nmax;
  double NCmin, NCmax;
  double Q[TYPE_COUNT][TYPE_COUNT];
  double alpha[TYPE_COUNT][TYPE_COUNT];
  double A[2][2];
  double BIJc_CC1, BIJc_CC2, BIJc_CC3;
  double BIJc_CH1, BIJc_CH2, BIJc_CH3;
  double BIJc_HH1, BIJc_HH2, BIJc_HH3;
  double BIJc[TYPE_COUNT][TYPE_COUNT][3];
  double Beta_CC1, Beta_CC2, Beta_CC3;
  double Beta_CH1, Beta_CH2, Beta_CH3;
  double Beta_HH1, Beta_HH2, Beta_HH3;
  double Beta[TYPE_COUNT][TYPE_COUNT][3];
  double rho[TYPE_COUNT][TYPE_COUNT];

  // b) czlon LJ
  double rcLJmin[TYPE_COUNT][TYPE_COUNT];
  double rcLJmax[TYPE_COUNT][TYPE_COUNT];
  double rcLJmaxsq[TYPE_COUNT][TYPE_COUNT];
  double bLJmin[TYPE_COUNT][TYPE_COUNT];
  double bLJmax[TYPE_COUNT][TYPE_COUNT];
  double epsilon[TYPE_COUNT][TYPE_COUNT];
  double sigma[TYPE_COUNT][TYPE_COUNT];

  // c) czlon torsyjny
  double thmin, thmax;
  double epsilonT_CCCC, epsilonT_CCCH, epsilonT_HCCH;
  double epsilonT[2][2];

  // d) splajny
  double gCdom[5], gC1[4][6], gC2[4][6];
  double gHdom[4], gH[3][6];
  double pCCdom[2][2], pCC[4][4][16];
  double pCHdom[2][2], pCH[4][4][16];
  double piCCdom[3][2], piCC[4][4][9][64];
  double piCHdom[3][2], piCH[4][4][9][64];
  double piHHdom[3][2], piHH[4][4][9][64];
  double Tijdom[3][2], Tijc[4][4][9][64];

  // splajny (polozenia wezlow)
  double PCCf[5][5], PCCdfdx[5][5], PCCdfdy[5][5];
  double PCHf[5][5], PCHdfdx[5][5], PCHdfdy[5][5];
  double piCCf[5][5][11], piCCdfdx[5][5][11], piCCdfdy[5][5][11], piCCdfdz[5][5][11];
  double piCHf[5][5][11], piCHdfdx[5][5][11], piCHdfdy[5][5][11], piCHdfdz[5][5][11];
  double piHHf[5][5][11], piHHdfdx[5][5][11], piHHdfdy[5][5][11], piHHdfdz[5][5][11];
  double Tf[5][5][10], Tdfdx[5][5][10], Tdfdy[5][5][10], Tdfdz[5][5][10];

  void readParameters( string file_name );
  void allocateMemory( );
  void deallocateMemory( );
  void initialize_constants( );
  void initialize_splines( );

  // ewaluacja splajnow
  double gSpline( double costh, double Nij, int typei );
  double PijSpline( double NijC, double NijH, int typei, int typej );
  double piRCSpline( double Nij, double Nji, double Nijconj, int typei, int typej );
  double TijSpline( double Nij, double Nji, double Nijconj );

  double Sp5th( double x, double *coeffs );
  double Spbicubic( double x, double y, double *coeffs );
  double Sptricubic( double x, double y, double z, double *coeffs );

  double bondorder( int i, int j, double *rji_vec, double rji );
  double bondorderLJ( int i, int j, double *rji_vec, double rji, double rji0 );

  void REBO_neighbours( );
  void E_REBO( );
  void E_LJ( );
  void E_TORSION( );

  void getLine( ifstream &file ) {
    string line;

    getline( file, line );
    if ( file.eof( ) ) {
      cout << "error in AIREBOForceField::(): end of file reached!" << endl;
      exit( 0 );
    }
  }

  void allocateREBO( int number_of_atoms ) {
    natoms = number_of_atoms;
    nC = new double[natoms]();
    nH = new double[natoms]();
    REBO_neighbours_num = new int[natoms];
    REBO_neighbours_list = new int*[natoms];
    REBO_neighbours_bonds = new vec3d*[natoms];
    for( int i = 0; i < natoms; i++ ) {
      REBO_neighbours_list[i] = new int[max_number_of_REBO_neighbours];
      REBO_neighbours_bonds[i] = new vec3d[max_number_of_REBO_neighbours];
    }
  }

  void deallocateREBO( ) {
    delete [] nC;
    delete [] nH;
    for( int i = 0; i < natoms; i++ ) {
      delete [] REBO_neighbours_list[i];
      delete [] REBO_neighbours_bonds[i];
    }
    delete [] REBO_neighbours_num;
    delete [] REBO_neighbours_list;
    delete [] REBO_neighbours_bonds;
  }

  double getLineAndConvertToDouble( ifstream &file ) {
    string line;
    istringstream iss;
    double value;

    getline( file, line );
    if ( file.eof( ) ) {
      cout << "error in AIREBOForceField::getLineAndConvertToDouble(): end of file reached!" << endl;
      exit( 0 );
    }

    iss.str( line );
    iss >> value;
    if ( value != value ) {
      cout << "error in AIREBOForceField::getLineAndConvertToDouble(): conversion error!" << endl;
      exit( 0 );
    }

    return value;
  }

  int getLineAndConvertToInt( ifstream &file ) {
    string line;
    istringstream iss;
    int value;

    getline( file, line );
    if ( file.eof( ) ) {
      cout << "error in AIREBOForceField::getLineAndConvertToInt(): end of file reached!" << endl;
      exit( 0 );
    }

    iss.str( line );
    iss >> value;
    if ( value != value ) {
      cout << "error in AIREBOForceField::getLineAndConvertToInt(): conversion error!" << endl;
      exit( 0 );
    }

    return value;
  }

  double min( double val1, double val2 ) {
    return ( val1 < val2 ) ? val1 : val2;
  }

  double max( double val1, double val2 ) {
    return ( val1 > val2 ) ? val1 : val2;
  }

  double square( double arg ) {
    return arg * arg;
  }

  double cube( double arg ) {
    return arg * arg * arg;
  }

  double pow4( double arg ) {
    return arg * arg * arg * arg;
  }

  double pow5( double arg ) {
    return arg * arg * arg * arg * arg;
  }

  double kronecker( const int a, const int b ) const {
    return ( a == b ) ? 1.0 : 0.0;
  }

  double SpRC( double Xij, int type1, int type2 ) const {
    double Xmin = rcMin[type1][type2];
    if ( Xij > rcMax[type1][type2] )
      return 0.0;
    if ( Xij <= Xmin )
      return 1.0;
    return 0.5 + 0.5 * cos( ( Xij - Xmin ) * pi_div_delta_rc[type1][type2] );
  }

  double SpRCP( double Xij, int type1, int type2 ) const {
    double Xmin = rcMin[type1][type2];
    if ( Xij > rcMaxP[type1][type2] )
      return 0.0;
    if ( Xij <= Xmin )
      return 1.0;
    return 0.5 + 0.5 * cos( ( Xij - Xmin ) * pi_div_delta_rcp[type1][type2] );
  }

  double SpN( double Xij, int type1, int type2 ) const {
    double Xmin = Nmin[type1][type2];
    if ( Xij > Nmax[type1][type2] )
      return 0.0;
    if ( Xij <= Xmin )
      return 1.0;
    return 0.5 + 0.5 * cos( ( Xij - Xmin ) * pi_div_delta_N[type1][type2] );
  }

  double SpNC( double Xij, int type1, int type2 ) const {
    double Xmin = NCmin[type1][type2];
    if ( Xij > NCmax[type1][type2] )
      return 0.0;
    if ( Xij <= Xmin )
      return 1.0;
    return 0.5 + 0.5 * cos( ( Xij - Xmin ) * pi_div_delta_NC[type1][type2] );
  }

  double Sp( double Xij, double Xmin, double Xmax ) const {
    // funkcja odciecia Sp
    double cutoff;
    double t = ( Xij - Xmin ) / ( Xmax - Xmin );
    // !!!nieoptymalnie, konieczne odraczanie dzielenia
    if ( t <= 0.0 )
      cutoff = 1.0;
    else if ( t >= 1.0 )
      cutoff = 0.0;
    else
      cutoff = 0.5 * ( 1.0 + cos( t * PI ) );
    return cutoff;
  }

  double Sp2( double Xij, double Xmin, double Xmax ) const {
    // funkcja odciecia Sp2
    double cutoff;
    double t = ( Xij - Xmin ) / ( Xmax - Xmin );
    // !!!nieoptymalnie, konieczne odraczanie dzielenia
    if ( t <= 0.0 )
      cutoff = 1.0;
    else if ( t >= 1.0 )
      cutoff = 0.0;
    else
      cutoff = ( 1.0 - ( t * t * ( 3.0 - 2.0 * t ) ) );
    return cutoff;
  }
};

#endif
