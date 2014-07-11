// based on LAMMPS implementation of AIREBO force field

#include "OpenREBO.h"

namespace OpenREBO {

  double AIREBO::compute( ) {
    assert( atom_list != nullptr );

    //atom_list->mixed_type_flag = true; //DEBUG
    if ( atom_list->is_pure_C( ) ) {
      REBO_neighbours_C( );
      energy_REBO = E_REBO_C( );
    } else {
      REBO_neighbours_CH( );
      energy_REBO = E_REBO_CH( );
    }
#ifndef LEAN_REBO
    energy_LJ = ( LJ_flag ) ? E_LJ( ) : 0;
    energy_torsion = ( torsion_flag ) ? E_Torsion( ) : 0;
#endif

    total_energy = energy_REBO + energy_LJ + energy_torsion;

    return total_energy;
  }

  void AIREBO::readNList( string filename ) {
    if ( atom_list != nullptr )
      delete atom_list;
    if ( rebo_atom_list != nullptr )
      delete rebo_atom_list;
    atom_list = new NList( filename, max_REBO_neighbours );
  }

  void AIREBO::readParameters( const string& file_name ) {
    ifstream file;
    file.open( file_name.c_str( ) );
    assert( file.is_open( ) );

    string line;
    do { // skip comments at beginning of file
      getline( file, line );
    } while( line.length( ) != 0 && line.at( 0 ) == '#' );

    rcMin[0][0] = getLineToDouble( file );
    rcMin[0][1] = getLineToDouble( file );
    rcMin[1][0] = rcMin[0][1];
    rcMin[1][1] = getLineToDouble( file );

    rcMax[0][0] = getLineToDouble( file );
    rcMax[0][1] = getLineToDouble( file );
    rcMax[1][0] = rcMax[0][1];
    rcMax[1][1] = getLineToDouble( file );

    rcMaxP[0][0] = getLineToDouble( file );
    rcMaxP[0][1] = getLineToDouble( file );
    rcMaxP[1][0] = rcMaxP[0][1];
    rcMaxP[1][1] = getLineToDouble( file );

    rcMaxSq[0][0] = rcMax[0][0] * rcMax[0][0];
    rcMaxSq[0][1] = rcMax[0][1] * rcMax[0][1];
    rcMaxSq[1][0] = rcMaxSq[0][1];
    rcMaxSq[1][1] = rcMax[1][1] * rcMax[1][1];

    smin = getLineToDouble( file );
    Nmin = getLineToDouble( file );
    Nmax = getLineToDouble( file );
    NCmin = getLineToDouble( file );
    NCmax = getLineToDouble( file );

    // MOD delt
    pi_div_delta_RC[0][0] = PI / ( rcMax[0][0] - rcMin[0][0] );
    pi_div_delta_RC[0][1] = PI / ( rcMax[0][1] - rcMin[0][1] );
    pi_div_delta_RC[1][0] = pi_div_delta_RC[0][1];
    pi_div_delta_RC[1][1] = PI / ( rcMax[1][1] - rcMin[1][1] );

    pi_div_delta_RCP[0][0] = PI / ( rcMaxP[0][0] - rcMin[0][0] );
    pi_div_delta_RCP[0][1] = PI / ( rcMaxP[0][1] - rcMin[0][1] );
    pi_div_delta_RCP[1][0] = pi_div_delta_RC[0][1];
    pi_div_delta_RCP[1][1] = PI / ( rcMaxP[1][1] - rcMin[1][1] );

    rcMin_CC = rcMin[0][0];
    rcMax_CC = rcMax[0][0];
    rcMaxP_CC = rcMaxP[0][0];
    pi_div_delta_RC_CC = pi_div_delta_RC[0][0];
    pi_div_delta_RCP_CC = pi_div_delta_RCP[0][0];

    pi_div_delta_N = PI / ( Nmax - Nmin );
    pi_div_delta_NC = PI / ( NCmax - NCmin );

    Q[0][0] = getLineToDouble( file );
    Q[0][1] = getLineToDouble( file );
    Q[1][0] = Q[0][1];
    Q[1][1] = getLineToDouble( file );
    alpha[0][0] = getLineToDouble( file );
    alpha[0][1] = getLineToDouble( file );
    alpha[1][0] = alpha[0][1];
    alpha[1][1] = getLineToDouble( file );
    A[0][0] = getLineToDouble( file );
    A[0][1] = getLineToDouble( file );
    A[1][0] = A[0][1];
    A[1][1] = getLineToDouble( file );

    BIJc[0][0][0] = getLineToDouble( file );
    BIJc[0][0][1] = getLineToDouble( file );
    BIJc[0][0][2] = getLineToDouble( file );
    BIJc[0][1][0] = getLineToDouble( file );
    BIJc[0][1][1] = getLineToDouble( file );
    BIJc[0][1][2] = getLineToDouble( file );
    BIJc[1][0][0] = BIJc[0][1][0];
    BIJc[1][0][1] = BIJc[0][1][1];
    BIJc[1][0][2] = BIJc[0][1][2];
    BIJc[1][1][0] = getLineToDouble( file );
    BIJc[1][1][1] = getLineToDouble( file );
    BIJc[1][1][2] = getLineToDouble( file );

    beta[0][0][0] = getLineToDouble( file );
    beta[0][0][1] = getLineToDouble( file );
    beta[0][0][2] = getLineToDouble( file );
    beta[0][1][0] = getLineToDouble( file );
    beta[0][1][1] = getLineToDouble( file );
    beta[0][1][2] = getLineToDouble( file );
    beta[1][0][0] = beta[0][1][0];
    beta[1][0][1] = beta[0][1][1];
    beta[1][0][2] = beta[0][1][2];
    beta[1][1][0] = getLineToDouble( file );
    beta[1][1][1] = getLineToDouble( file );
    beta[1][1][2] = getLineToDouble( file );

    rho[0][0] = getLineToDouble( file );
    rho[0][1] = getLineToDouble( file );
    rho[1][0] = rho[0][1];
    rho[1][1] = getLineToDouble( file );

    // b) czlon LJ
    rcLJmin[0][0] = getLineToDouble( file );
    rcLJmin[0][1] = getLineToDouble( file );
    rcLJmin[1][0] = rcLJmin[0][1];
    rcLJmin[1][1] = getLineToDouble( file );

    rcLJmax[0][0] = getLineToDouble( file );
    rcLJmax[0][1] = getLineToDouble( file );
    rcLJmax[1][0] = rcLJmax[0][1];
    rcLJmax[1][1] = getLineToDouble( file );

    bLJmin[0][0] = getLineToDouble( file );
    bLJmin[0][1] = getLineToDouble( file );
    bLJmin[1][0] = bLJmin[0][1];
    bLJmin[1][1] = getLineToDouble( file );

    bLJmax[0][0] = getLineToDouble( file );
    bLJmax[0][1] = getLineToDouble( file );
    bLJmin[1][0] = bLJmin[0][1];
    bLJmax[1][1] = getLineToDouble( file );

    inv_delta_bLJ[0][0] = PI / ( bLJmax[0][0] - bLJmin[0][0] );
    inv_delta_bLJ[0][1] = PI / ( bLJmax[0][1] - bLJmin[0][1] );
    inv_delta_bLJ[1][0] = inv_delta_bLJ[0][1];
    inv_delta_bLJ[1][1] = PI / ( bLJmax[1][1] - bLJmin[1][1] );

    rcLJmaxsq[0][0] = rcLJmax[0][0] * rcLJmax[0][0];
    rcLJmaxsq[0][1] = rcLJmax[0][1] * rcLJmax[0][1];
    rcLJmaxsq[1][0] = rcLJmaxsq[0][1];
    rcLJmaxsq[1][1] = rcLJmax[1][1] * rcLJmax[1][1];

    epsilon[0][0] = getLineToDouble( file );
    epsilon[0][1] = getLineToDouble( file );
    epsilon[1][0] = epsilon[0][1];
    epsilon[1][1] = getLineToDouble( file );

    sigma[0][0] = getLineToDouble( file );
    sigma[0][1] = getLineToDouble( file );
    sigma[1][0] = sigma[0][1];
    sigma[1][1] = getLineToDouble( file );

    epsilonT[0][0] = getLineToDouble( file ); // CCCC
    epsilonT[0][1] = getLineToDouble( file ); // CCCH
    epsilonT[1][0] = epsilonT[0][1]; // CCCH
    epsilonT[1][1] = getLineToDouble( file ); // HCCH

    readSplineGC( file, gCdom, gC1, gC2 );
    readSplineGH( file, gHdom, gH );
    readSplineP( file, pCCdom, pCC );
    readSplineP( file, pCHdom, pCH );

    readSplinePI( file, piCCdom, piCC );
    readSplinePI( file, piCHdom, piCH );
    readSplinePI( file, piHHdom, piHH );
    readSplinePI( file, Tijdom, Tijc );

    file.close( );

    thmin = -1.0;
    thmax = -0.995;
    inv_th_delta = 1.0 / ( thmax - thmin );
  }

  void AIREBO::readSplineGC( ifstream& ifs, double Vdom[5], double v1[4][6], double v2[4][6] ) {
    for( int i = 0; i < 3; i++ )
      getLine( ifs );
    int nr_domains = getLineToInt( ifs );
    for( int i = 0; i < nr_domains; i++ )
      Vdom[i] = getLineToDouble( ifs );
    //assert(Vdom[i] > Vdom[i]-1);
    getLine( ifs );
    nr_domains--;
    for( int i = 0; i < nr_domains; i++ )
      for( int j = 0; j < 6; j++ )
        v1[i][j] = getLineToDouble( ifs );
    getLine( ifs );
    for( int i = 0; i < nr_domains; i++ )
      for( int j = 0; j < 6; j++ )
        v2[i][j] = getLineToDouble( ifs );
  }

  void AIREBO::readSplineGH( ifstream& ifs, double Vdom[4], double v[3][6] ) {
    for( int i = 0; i < 3; i++ )
      getLine( ifs );
    int nr_domains = getLineToInt( ifs );
    for( int i = 0; i < nr_domains; i++ )
      Vdom[i] = getLineToDouble( ifs );
    //assert(Vdom[i] > Vdom[i]-1);
    getLine( ifs );
    nr_domains--;
    for( int i = 0; i < nr_domains; i++ )
      for( int j = 0; j < 6; j++ )
        v[i][j] = getLineToDouble( ifs );
  }

  void AIREBO::readSplineP( ifstream& ifs, double Vdom[2][2], double v[4][4][16] ) {
    for( int i = 0; i < 3; i++ )
      getLine( ifs );
    int nr_domains = getLineToInt( ifs );
    int max = nr_domains / 2;
    for( int i = 0; i < max; i++ )
      for( int j = 0; j < max; j++ )
        Vdom[i][j] = getLineToDouble( ifs );
    //assert(Vdom[i][j] > Vdom[i][j-1]);
    getLine( ifs );
    for( int i = 0, max_i = (int) Vdom[0][1]; i < max_i; i++ )
      for( int j = 0, max_j = (int) Vdom[1][1]; j < max_j; j++ )
        for( int k = 0; k < 16; k++ )
          v[i][j][k] = getLineToDouble( ifs );
  }

  void AIREBO::readSplinePI( ifstream& ifs, double Vdom[3][2], double v[4][4][9][64] ) {
    for( int i = 0; i < 3; i++ )
      getLine( ifs );
    int nr_domains = getLineToInt( ifs );
    for( int i = 0, max_i = nr_domains / 2; i < max_i; i++ )
      for( int j = 0, max_j = nr_domains / 3; j < max_j; j++ )
        Vdom[i][j] = getLineToDouble( ifs );
    //assert(Vdom[i][j] > Vdom[i][j-1]);
    getLine( ifs );
    for( int i = 0, max_i = (int) Vdom[0][1]; i < max_i; i++ )
      for( int j = 0, max_j = (int) Vdom[1][1]; j < max_j; j++ )
        for( int k = 0, max_k = (int) Vdom[2][1]; k < max_k; k++ )
          for( int l = 0; l < 64; l++ )
            v[i][j][k][l] = getLineToDouble( ifs );
  }

  void AIREBO::initializeSplines( ) {
    int i, j, k;

    for( i = 0; i < 5; i++ )
      for( j = 0; j < 5; j++ ) {
        PCCf[i][j] = 0.0;
        PCCdfdx[i][j] = 0.0;
        PCCdfdy[i][j] = 0.0;
        PCHf[i][j] = 0.0;
        PCHdfdx[i][j] = 0.0;
        PCHdfdy[i][j] = 0.0;
      }

    PCCf[0][2] = -0.00050;
    PCCf[0][3] = 0.0161253646;
    PCCf[1][1] = -0.010960;
    PCCf[1][2] = 0.00632624824;
    PCCf[2][0] = -0.0276030;
    PCCf[2][1] = 0.00317953083;

    PCHf[0][1] = 0.209336733;
    PCHf[0][2] = -0.0644496154;
    PCHf[0][3] = -0.303927546;
    PCHf[1][0] = 0.010;
    PCHf[1][1] = -0.125123401;
    PCHf[1][2] = -0.298905246;
    PCHf[2][0] = -0.122042146;
    PCHf[2][1] = -0.300529172;
    PCHf[3][0] = -0.307584705;

    for( i = 0; i < 5; i++ )
      for( j = 0; j < 5; j++ )
        for( k = 0; k < 10; k++ ) {
          piCCf[i][j][k] = 0.0;
          piCCdfdx[i][j][k] = 0.0;
          piCCdfdy[i][j][k] = 0.0;
          piCCdfdz[i][j][k] = 0.0;
          piCHf[i][j][k] = 0.0;
          piCHdfdx[i][j][k] = 0.0;
          piCHdfdy[i][j][k] = 0.0;
          piCHdfdz[i][j][k] = 0.0;
          piHHf[i][j][k] = 0.0;
          piHHdfdx[i][j][k] = 0.0;
          piHHdfdy[i][j][k] = 0.0;
          piHHdfdz[i][j][k] = 0.0;
          Tf[i][j][k] = 0.0;
          Tdfdx[i][j][k] = 0.0;
          Tdfdy[i][j][k] = 0.0;
          Tdfdz[i][j][k] = 0.0;
        }

    for( i = 3; i < 10; i++ )
      piCCf[0][0][i] = 0.0049586079;
    piCCf[1][0][1] = 0.021693495;
    piCCf[0][1][1] = 0.021693495;
    for( i = 2; i < 10; i++ )
      piCCf[1][0][i] = 0.0049586079;
    for( i = 2; i < 10; i++ )
      piCCf[0][1][i] = 0.0049586079;
    piCCf[1][1][1] = 0.05250;
    piCCf[1][1][2] = -0.002088750;
    for( i = 3; i < 10; i++ )
      piCCf[1][1][i] = -0.00804280;
    piCCf[2][0][1] = 0.024698831850;
    piCCf[0][2][1] = 0.024698831850;
    piCCf[2][0][2] = -0.00597133450;
    piCCf[0][2][2] = -0.00597133450;
    for( i = 3; i < 10; i++ )
      piCCf[2][0][i] = 0.0049586079;
    for( i = 3; i < 10; i++ )
      piCCf[0][2][i] = 0.0049586079;
    piCCf[2][1][1] = 0.00482478490;
    piCCf[1][2][1] = 0.00482478490;
    piCCf[2][1][2] = 0.0150;
    piCCf[1][2][2] = 0.0150;
    piCCf[2][1][3] = -0.010;
    piCCf[1][2][3] = -0.010;
    piCCf[2][1][4] = -0.01168893870;
    piCCf[1][2][4] = -0.01168893870;
    piCCf[2][1][5] = -0.013377877400;
    piCCf[1][2][5] = -0.013377877400;
    piCCf[2][1][6] = -0.015066816000;
    piCCf[1][2][6] = -0.015066816000;
    for( i = 7; i < 10; i++ )
      piCCf[2][1][i] = -0.015066816000;
    for( i = 7; i < 10; i++ )
      piCCf[1][2][i] = -0.015066816000;
    piCCf[2][2][1] = 0.0472247850;
    piCCf[2][2][2] = 0.0110;
    piCCf[2][2][3] = 0.0198529350;
    piCCf[2][2][4] = 0.01654411250;
    piCCf[2][2][5] = 0.013235290;
    piCCf[2][2][6] = 0.00992646749999;
    piCCf[2][2][7] = 0.006617644999;
    piCCf[2][2][8] = 0.00330882250;
    piCCf[3][0][1] = -0.05989946750;
    piCCf[0][3][1] = -0.05989946750;
    piCCf[3][0][2] = -0.05989946750;
    piCCf[0][3][2] = -0.05989946750;
    for( i = 3; i < 10; i++ )
      piCCf[3][0][i] = 0.0049586079;
    for( i = 3; i < 10; i++ )
      piCCf[0][3][i] = 0.0049586079;
    piCCf[3][1][2] = -0.0624183760;
    piCCf[1][3][2] = -0.0624183760;
    for( i = 3; i < 10; i++ )
      piCCf[3][1][i] = -0.0624183760;
    for( i = 3; i < 10; i++ )
      piCCf[1][3][i] = -0.0624183760;
    piCCf[3][2][1] = -0.02235469150;
    piCCf[2][3][1] = -0.02235469150;
    for( i = 2; i < 10; i++ )
      piCCf[3][2][i] = -0.02235469150;
    for( i = 2; i < 10; i++ )
      piCCf[2][3][i] = -0.02235469150;

    piCCdfdx[2][1][1] = -0.026250;
    piCCdfdx[2][1][5] = -0.0271880;
    piCCdfdx[2][1][6] = -0.0271880;
    for( i = 7; i < 10; i++ )
      piCCdfdx[2][1][i] = -0.0271880;
    piCCdfdx[1][3][2] = 0.0187723882;
    for( i = 2; i < 10; i++ )
      piCCdfdx[2][3][i] = 0.031209;

    piCCdfdy[1][2][1] = -0.026250;
    piCCdfdy[1][2][5] = -0.0271880;
    piCCdfdy[1][2][6] = -0.0271880;
    for( i = 7; i < 10; i++ )
      piCCdfdy[1][2][i] = -0.0271880;
    piCCdfdy[3][1][2] = 0.0187723882;
    for( i = 2; i < 10; i++ )
      piCCdfdy[3][2][i] = 0.031209;

    piCCdfdz[1][1][2] = -0.0302715;
    piCCdfdz[2][1][4] = -0.0100220;
    piCCdfdz[1][2][4] = -0.0100220;
    piCCdfdz[2][1][5] = -0.0100220;
    piCCdfdz[1][2][5] = -0.0100220;
    for( i = 4; i < 9; i++ )
      piCCdfdz[2][2][i] = -0.0033090;

    // make top end of piCC flat instead of zero also enforces some symmetry
    i = 4;
    for( j = 0; j < 4; j++ )
      for( k = 1; k < 11; k++ )
        piCCf[i][j][k] = piCCf[i - 1][j][k];
    for( i = 0; i < 4; i++ )
      for( j = i + 1; j < 5; j++ )
        for( k = 1; k < 11; k++ )
          piCCf[i][j][k] = piCCf[j][i][k];
    for( k = 1; k < 11; k++ )
      piCCf[4][4][k] = piCCf[3][4][k];
    k = 10;
    for( i = 0; i < 5; i++ )
      for( j = 0; j < 5; j++ )
        piCCf[i][j][k] = piCCf[i][j][k - 1];

    piCHf[1][1][1] = -0.050;
    piCHf[1][1][2] = -0.050;
    piCHf[1][1][3] = -0.30;
    for( i = 4; i < 10; i++ )
      piCHf[1][1][i] = -0.050;
    for( i = 5; i < 10; i++ )
      piCHf[2][0][i] = -0.004523893758064;
    for( i = 5; i < 10; i++ )
      piCHf[0][2][i] = -0.004523893758064;
    piCHf[2][1][2] = -0.250;
    piCHf[1][2][2] = -0.250;
    piCHf[2][1][3] = -0.250;
    piCHf[1][2][3] = -0.250;
    piCHf[3][1][1] = -0.10;
    piCHf[1][3][1] = -0.10;
    piCHf[3][1][2] = -0.125;
    piCHf[1][3][2] = -0.125;
    piCHf[3][1][3] = -0.125;
    piCHf[1][3][3] = -0.125;
    for( i = 4; i < 10; i++ )
      piCHf[3][1][i] = -0.10;
    for( i = 4; i < 10; i++ )
      piCHf[1][3][i] = -0.10;

    // make top end of piCH flat instead of zero also enforces some symmetry
    i = 4;
    for( j = 0; j < 4; j++ )
      for( k = 1; k < 11; k++ )
        piCHf[i][j][k] = piCHf[i - 1][j][k];
    for( i = 0; i < 4; i++ )
      for( j = i + 1; j < 5; j++ )
        for( k = 1; k < 11; k++ )
          piCHf[i][j][k] = piCHf[j][i][k];
    for( k = 1; k < 11; k++ )
      piCHf[4][4][k] = piCHf[3][4][k];
    k = 10;
    for( i = 0; i < 5; i++ )
      for( j = 0; j < 5; j++ )
        piCHf[i][j][k] = piCHf[i][j][k - 1];

    piHHf[1][1][1] = 0.124915958;

    Tf[2][2][1] = -0.035140;
    for( i = 2; i < 10; i++ )
      Tf[2][2][i] = -0.0040480;

    initializeSplines2( );
  }

  void AIREBO::initializeSplines2( ) {
    gSplineC1_low = Sp5th( gCdom[0], gC1[0] );
    gSplineC1_hi = Sp5th( gCdom[4], gC1[3] );
    gSplineC2_low = Sp5th( gCdom[0], gC2[0] );
    gSplineC2_hi = Sp5th( gCdom[4], gC2[3] );
    gSplineH_low = Sp5th( gHdom[0], gH[0] );
    gSplineH_hi = Sp5th( gHdom[3], gH[2] );
    gSplineC_low_delta = gSplineC1_low - gSplineC2_low;
    gSplineC_hi_delta = gSplineC1_hi - gSplineC2_hi;
  }

}