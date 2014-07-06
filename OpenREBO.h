/*
    File:   OpenREBO.h
    Author: vaxquis
    based on LAMMPS implementation of AIREBO potential
 */
#pragma once

#ifndef AIREBO_FORCE_FIELD_H
#define AIREBO_FORCE_FIELD_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cassert>

#define LEAN_REBO

#include "nlist.h"

namespace OpenREBO {
  using std::ifstream;
  using std::getline;
  using std::string;
  using std::endl;
  using std::istringstream;
  using std::cout;

  const int SPL_DIM = 4;

  class AIREBO {
    friend class NList; // for the sake of SpXX

  public:

    AIREBO( const string& filename, double cutR_LJ_sigma,
            bool LJ_flag, bool torsion_flag, int max_REBO_neighbours ) {
      atom_list = nullptr;
      rebo_atom_list = nullptr;

      total_energy = 0.0;
      energy_REBO = 0.0;
      energy_LJ = 0.0;
      energy_torsion = 0.0;

      readParameters( filename );
      this->LJ_flag = LJ_flag;
      if ( LJ_flag ) {
        assert( cutR_LJ_sigma > 0.0 );
        calcCutMax( cutR_LJ_sigma );
      }
      this->torsion_flag = torsion_flag;

      assert( max_REBO_neighbours >= 1 );
      this->max_REBO_neighbours = max_REBO_neighbours;

      initializeSplines( );
    }

    virtual ~AIREBO( ) {
      if ( atom_list != nullptr )
        delete atom_list;
      if ( rebo_atom_list != nullptr )
        delete rebo_atom_list;
    };

    double getCutoffRadius( ) {
      return cutMax;
    }

    double getEnergyTotal( ) {
      return total_energy;
    }

    double getEnergyREBO( ) {
      return energy_REBO;
    }

    double getEnergyLJ( ) {
      return energy_LJ;
    }

    double getEnergyTorsion( ) {
      return energy_torsion;
    }

    double compute( );
    void readNList( string filename );


  private:
    bool LJ_flag, torsion_flag;
    int max_REBO_neighbours;

    double total_energy, energy_REBO, energy_LJ, energy_torsion;

    NList* atom_list;
    RNList* rebo_atom_list;

    double cutLJsq[2][2];
    double lj1[2][2], lj2[2][2], lj3[2][2], lj4[2][2];
    double cut3rebo, cutLJrebo, cutLJrebosq;
    double cutMax;

    // basic REBO vars
    double rcMin[TYPE_COUNT][TYPE_COUNT], rcMax[TYPE_COUNT][TYPE_COUNT];
    double rcMaxSq[TYPE_COUNT][TYPE_COUNT], rcMaxP[TYPE_COUNT][TYPE_COUNT];
    double pi_div_delta_RC[TYPE_COUNT][TYPE_COUNT], pi_div_delta_RCP[TYPE_COUNT][TYPE_COUNT];

    double smin;
    double Nmin, Nmax;
    double NCmin, NCmax;
    double pi_div_delta_N, pi_div_delta_NC;
    double Q[TYPE_COUNT][TYPE_COUNT];
    double A[2][2];
    double alpha[TYPE_COUNT][TYPE_COUNT];
    double beta[TYPE_COUNT][TYPE_COUNT][3];
    double BIJc[TYPE_COUNT][TYPE_COUNT][3];
    double rho[TYPE_COUNT][TYPE_COUNT];

    // LJ vars
    double rcLJmin[TYPE_COUNT][TYPE_COUNT], rcLJmax[TYPE_COUNT][TYPE_COUNT];
    double rcLJmaxsq[TYPE_COUNT][TYPE_COUNT];
    double bLJmin[TYPE_COUNT][TYPE_COUNT], bLJmax[TYPE_COUNT][TYPE_COUNT];
    double inv_delta_bLJ[TYPE_COUNT][TYPE_COUNT];
    double epsilon[TYPE_COUNT][TYPE_COUNT];
    double sigma[TYPE_COUNT][TYPE_COUNT];

    // torsion vars
    double thmin, thmax, inv_th_delta;
    double epsilonT[TYPE_COUNT][TYPE_COUNT];

    // splines
    double gCdom[5], gC1[4][6], gC2[4][6];
    double gHdom[4], gH[3][6];
    double pCCdom[2][2], pCC[4][4][16];
    double pCHdom[2][2], pCH[4][4][16];
    double piCCdom[3][2], piCC[4][4][9][64];
    double piCHdom[3][2], piCH[4][4][9][64];
    double piHHdom[3][2], piHH[4][4][9][64];
    double Tijdom[3][2], Tijc[4][4][9][64];

    // spline nodes
    double PCCf[5][5], PCCdfdx[5][5], PCCdfdy[5][5];
    double PCHf[5][5], PCHdfdx[5][5], PCHdfdy[5][5];
    double piCCf[5][5][11], piCCdfdx[5][5][11], piCCdfdy[5][5][11], piCCdfdz[5][5][11];
    double piCHf[5][5][11], piCHdfdx[5][5][11], piCHdfdy[5][5][11], piCHdfdz[5][5][11];
    double piHHf[5][5][11], piHHdfdx[5][5][11], piHHdfdy[5][5][11], piHHdfdz[5][5][11];
    double Tf[5][5][11], Tdfdx[5][5][11], Tdfdy[5][5][11], Tdfdz[5][5][11]; // 3rd dim is 10 actually, this is only for compiler ease

    // spline helpers
    double gSplineC1_low, gSplineC1_hi;
    double gSplineC2_low, gSplineC2_hi;
    double gSplineH_low, gSplineH_hi;
    double gSplineC_low_delta, gSplineC_hi_delta;

    void readSplineGC( ifstream& ifs, double Vdom[5], double v1[4][6], double v2[4][6] );
    void readSplineGH( ifstream& ifs, double Vdom[4], double v[3][6] );
    void readSplineP( ifstream& ifs, double Vdom[2][2], double v[4][4][16] );
    void readSplinePI( ifstream& ifs, double Vdom[3][2], double v[4][4][9][64] );

    void readParameters( const string& file_name );
    void initializeSplines( );
    void initializeSplines2( );

    void allocateMemory( );
    void deallocateMemory( );

    double gSpline( double costh, double Nij, int typei );
    double gSplineC( double costh, double Nij );
    double gSplineH( double costh );

    double PijSpline( double NijC, double NijH, int typei, int typej );
    double PijSplineC_only( double NijC );

    double TricubicSplineMain( double Nij, double Nji, double Nijconj,
                               double dom[3][2], double coeff_f[5][5][11], double coeff[SPL_DIM][SPL_DIM][9][SPL_DIM*SPL_DIM*SPL_DIM] );

    double piRCSpline( double Nij, double Nji, double Nijconj, int typei, int typej );
    double piRCSplineCC( double Nij, double Nji, double Nijconj );
    double piRCSplineCH( double Nij, double Nji, double Nijconj );
    double piRCSplineHH( double Nij, double Nji, double Nijconj );

    double TijSpline( double Nij, double Nji, double Nijconj );

    double Sp5th( double x, double *coeffs );
    double Spbicubic( double x, double y, double *coeffs );
    double Sptricubic( double x, double y, double z, double *coeffs );

    void REBO_neighbours_C( );
    double E_REBO_C( );
    double bond_order_C( int i, int j, const bond* b );
    double bond_order_1_C( int a1, int a2, const double r12_vec[3], double r12, double inv_r12,
                           double NC, double &p, const neighbor_tracker *nt );

    void REBO_neighbours_CH( );
    double E_REBO_CH( );
    double bond_order_CH( int i, int j, const bond* b );
    double bond_order_1_CH( int a1, int a2, int type1, int type2, const double r12_vec[3], double r12, double inv_r12,
                            double NC, double NH, double NTotal, double &p, const neighbor_tracker *nt );

#ifndef LEAN_REBO
    double E_LJ( );
    double E_Torsion( );
    double bond_orderLJ( int i, int j, double *rji_vec, double rji, double rji0 );
#endif

    string getLine( ifstream& ifs ) {
      string line;
      getline( ifs, line );
      assert( !ifs.eof( ) );
      return line;
    }

    double getLineToDouble( ifstream& ifs ) {
      istringstream iss( getLine( ifs ) );
      double value;
      iss >> value;
      assert( value == value ); // !NaN
      return value;
    }

    int getLineToInt( ifstream& ifs ) {
      istringstream iss( getLine( ifs ) );
      int value;
      iss >> value;
      return value;
    }

    void calcLJpair( int i, int j, double cutR_LJ_sigmas ) {
      double tmp = pow( sigma[i][j], 6.0 );
      lj4[i][j] = 4.0 * epsilon[i][j] * tmp;
      lj3[i][j] = lj4[i][j] * tmp;
      lj2[i][j] = lj4[i][j] * 6.0;
      lj1[i][j] = lj3[i][j] * 12.0;
      tmp = cutR_LJ_sigmas * sigma[i][j];
      cutLJsq[i][j] = tmp * tmp;
    }

    void calcCutMax( double cutR_LJ_sigmas ) {
      calcLJpair( 0, 0, cutR_LJ_sigmas );
      calcLJpair( 0, 1, cutR_LJ_sigmas );
      calcLJpair( 1, 1, cutR_LJ_sigmas );
      lj4[1][0] = lj4[0][1];
      lj3[1][0] = lj3[0][1];
      lj2[1][0] = lj2[0][1];
      lj1[1][0] = lj1[0][1];
      cutLJsq[1][0] = cutLJsq[0][1];

      cut3rebo = 3.0 * rcMax[0][0];
      cutLJrebo = rcLJmax[0][0] + rcMax[0][0];
      cutLJrebosq = cutLJrebo * cutLJrebo;

      cutMax = cut3rebo;

      if ( !LJ_flag )
        return;

      double tmp_cutmax = rcLJmax[0][0] + 2.0 * rcMax[0][0];
      if ( tmp_cutmax > cutMax )
        cutMax = tmp_cutmax;
      tmp_cutmax = cutR_LJ_sigmas * sigma[0][0];
      if ( tmp_cutmax > cutMax )
        cutMax = tmp_cutmax;
    }

    double SpRC( double Xij, int type1, int type2 ) const {
      double Xmin = rcMin[type1][type2];
      if ( Xij >= rcMax[type1][type2] )
        return 0.0;
      if ( Xij <= Xmin )
        return 1.0;
      return 0.5 + 0.5 * cos( ( Xij - Xmin ) * pi_div_delta_RC[type1][type2] );
    }

    double SpRCP( double Xij, int type1, int type2 ) const {
      double Xmin = rcMin[type1][type2];
      if ( Xij >= rcMaxP[type1][type2] )
        return 0.0;
      if ( Xij <= Xmin )
        return 1.0;
      return 0.5 + 0.5 * cos( ( Xij - Xmin ) * pi_div_delta_RCP[type1][type2] );
    }

    double SpN( double Xij ) const {
      if ( Xij >= Nmax )
        return 0.0;
      if ( Xij <= Nmin )
        return 1.0;
      return 0.5 + 0.5 * cos( ( Xij - Nmin ) * pi_div_delta_N );
    }

    double SpNC( double Xij ) const {
      if ( Xij >= NCmax )
        return 0.0;
      if ( Xij <= NCmin )
        return 1.0;
      return 0.5 + 0.5 * cos( ( Xij - NCmin ) * pi_div_delta_NC );
    }

    double Sp2th( double Xij ) const {
      if ( Xij <= thmin )
        return 1.0;
      if ( Xij >= thmax )
        return 0.0;

      double t = ( Xij - thmin ) * inv_th_delta;

      return 1.0 - t * t * ( 3.0 - 2.0 * t );
    }

    double Sp2bLJ( double Xij, int type1, int type2 ) const {
      double Xmin = bLJmin[type1][type2];
      if ( Xij >= bLJmax[type1][type2] )
        return 0.0;
      if ( Xij <= Xmin )
        return 1.0;

      double t = ( Xij - Xmin ) * inv_delta_bLJ[type1][type2];

      return 1.0 - t * t * ( 3.0 - 2.0 * t );
    }
  };

}

#endif /* AIREBO_FORCE_FIELD_H */
