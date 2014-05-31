// written by Szymon Winczewski
// based on LAMMPS implementation of AIREBO force field

#include "airebo_force_field.h"

using namespace std;

AIREBOForceField::AIREBOForceField( string file_name, double cutlj,
        bool ljflag, bool torflag, int max_number_of_REBO_neighbours ) {
  natoms = 0;
  type = NULL;
  nC = NULL;
  nH = NULL;
  cutljsq = NULL;
  lj1 = NULL;
  lj2 = NULL;
  lj3 = NULL;
  lj4 = NULL;

  total_energy = 0.0;
  energy_rebo = 0.0;
  energy_lj = 0.0;
  energy_torsion = 0.0;

  readParameters( file_name );
  if ( cutlj <= 0.0 ) {
    cout << "error in AIREBOForceField::AIREBOForceField(): incorrect parameter (cutlj)!" << endl;
    exit( 0 );
  }
  this->cutlj = cutlj;
  this->ljflag = ljflag;
  this->torflag = torflag;

  if ( max_number_of_REBO_neighbours < 1 ) {
    cout << "error in AIREBOForceField::AIREBOForceField(): incorrect parameter (max_number_of_REBO_neighbours)!" << endl;
    exit( 0 );
  }
  this->max_number_of_REBO_neighbours = max_number_of_REBO_neighbours;

  allocateMemory( );
  initialize_constants( );
  initialize_splines( );
}

AIREBOForceField::~AIREBOForceField( ) {
  deallocateMemory( );
}

double AIREBOForceField::getCutoffRadius( ) {
  return cutmax;
}

double AIREBOForceField::compute( int number_of_atoms, int *type,
        int *neighbours_num, int **neighbours_list,
        vec3d **neighbours_bonds ) {
#ifdef AIREBO_DEBUG
  cout << "AIREBOForceField::compute() started!" << endl;
  cout << "   para_number_of_atoms    " << number_of_atoms << endl;
  cout << "   para_type               " << type << endl;
  cout << "   para_neighbours_num     " << neighbours_num << endl;
  cout << "   para_neighbours_list    " << neighbours_list << endl;
  cout << "   para_neighbours_bonds   " << neighbours_bonds << endl;
  cout << endl;
#endif

  if ( number_of_atoms <= 0 ) {
    cout << "error in AIREBOForceField::compute(): incorrect parameter (number_of_atoms)" << endl;
    exit( 0 );
  }

  if ( natoms == 0 )
    allocateREBO( number_of_atoms );
  else if ( natoms != number_of_atoms ) {
    deallocateREBO( );
    allocateREBO( number_of_atoms );
  }

  if ( type == NULL ) {
    cout << "error in AIREBOForceField::compute(): incorrect (NULL) parameter (type)" << endl;
    exit( 0 );
  }
  this->type = type;

  if ( neighbours_num == NULL ) {
    cout << "error in AIREBOForceField::compute(): incorrect (NULL) parameter (para_neighbours_num)" << endl;
    exit( 0 );
  }
  this->neighbours_num = neighbours_num;

  if ( neighbours_list == NULL ) {
    cout << "error in AIREBOForceField::compute(): incorrect (NULL) parameter (para_neighbours_list)" << endl;
    exit( 0 );
  }
  this->neighbours_list = neighbours_list;

  if ( neighbours_bonds == NULL ) {
    cout << "error in AIREBOForceField::compute(): incorrect (NULL) parameter (para_neighbours_bonds)" << endl;
    exit( 0 );
  }
  this->neighbours_bonds = neighbours_bonds;

  total_energy = 0.0;
  energy_rebo = 0.0;
  energy_lj = 0.0;
  energy_torsion = 0.0;

#ifdef AIREBO_DEBUG
  cout << "starting REBO_neighbours()!" << endl;
#endif

  REBO_neighbours( );

#ifdef AIREBO_DEBUG
  cout << "REBO_neighbours() finished!" << endl;
#endif

#ifdef AIREBO_DEBUG
  cout << "starting E_REBO()!" << endl;
#endif

  E_REBO( );

#ifdef AIREBO_DEBUG
  cout << "E_REBO() finished!" << endl;
  cout << "energy_rebo = " << energy_rebo << endl;
  cout << endl;
#endif

  if ( ljflag ) {
#ifdef AIREBO_DEBUG
    cout << "starting E_LJ()!" << endl;
#endif

    E_LJ( );

#ifdef AIREBO_DEBUG
    cout << "E_LJ() finished!" << endl;
    cout << "energy_lj = " << energy_lj << endl;
    cout << endl;
#endif
  }

  if ( torflag ) {
#ifdef AIREBO_DEBUG
    cout << "starting E_TORSION()!" << endl;
#endif

    E_TORSION( );

#ifdef AIREBO_DEBUG
    cout << "E_TORSION() finished!" << endl;
    cout << "energy_torsion = " << energy_torsion << endl;
    cout << endl;
#endif
  }

  total_energy = energy_rebo + energy_lj + energy_torsion;

  return total_energy;
}

double AIREBOForceField::getTotalEnergy( ) {
  return total_energy;
}

double AIREBOForceField::getREBOEnergy( ) {
  return energy_rebo;
}

double AIREBOForceField::getLJEnergy( ) {
  return energy_lj;
}

double AIREBOForceField::getTORSIONEnergy( ) {
  return energy_torsion;
}

void AIREBOForceField::readParameters( string file_name ) {
  ifstream file;
  string line;
  int i, j, k, l;
  int number_of_domains;

  file.open( file_name.c_str( ) );
  if ( !file.is_open( ) ) {
    cout << "error in AIREBOForceField::readParameters(): could not locate AIREBO potential file '" << file_name << "'!" << endl;
    exit( 0 );
  }

  do { // skip comments at beginning of file
    getline( file, line );
  } while( line.length( ) != 0 && line.at( 0 ) == '#' );

  // a) czlon REBO
  rcMin[0][0] = getLineAndConvertToDouble( file );
  rcMin[0][1] = getLineAndConvertToDouble( file );
  rcMin[1][0] = rcMin[0][1];
  rcMin[1][1] = getLineAndConvertToDouble( file );

  rcMax[0][0] = getLineAndConvertToDouble( file );
  rcMax[0][1] = getLineAndConvertToDouble( file );
  rcMax[1][0] = rcMax[0][1];
  rcMax[1][1] = getLineAndConvertToDouble( file );

  rcMaxP[0][0] = getLineAndConvertToDouble( file );
  rcMaxP[0][1] = getLineAndConvertToDouble( file );
  rcMaxP[1][0] = rcMaxP[0][1];
  rcMaxP[1][1] = getLineAndConvertToDouble( file );

  rcMaxSq[0][0] = rcMax[0][0] * rcMax[0][0];
  rcMaxSq[0][1] = rcMax[0][1] * rcMax[0][1];
  rcMaxSq[1][0] = rcMaxSq[0][1];
  rcMaxSq[1][1] = rcMax[1][1] * rcMax[1][1];

  smin = getLineAndConvertToDouble( file );
  Nmin = getLineAndConvertToDouble( file );
  Nmax = getLineAndConvertToDouble( file );
  NCmin = getLineAndConvertToDouble( file );
  NCmax = getLineAndConvertToDouble( file );

  // MOD delt
  pi_div_delta_rc[0][0] = PI / ( rcMax[0][0] - rcMin[0][0] );
  pi_div_delta_rc[0][1] = PI / ( rcMax[0][1] - rcMin[0][1] );
  pi_div_delta_rc[1][0] = pi_div_delta_rc[0][1];
  pi_div_delta_rc[1][1] = PI / ( rcMax[1][1] - rcMin[1][1] );

  pi_div_delta_rcp[0][0] = PI / ( rcMaxP[0][0] - rcMin[0][0] );
  pi_div_delta_rcp[0][1] = PI / ( rcMaxP[0][1] - rcMin[0][1] );
  pi_div_delta_rcp[1][0] = pi_div_delta_rc[0][1];
  pi_div_delta_rcp[1][1] = PI / ( rcMaxP[1][1] - rcMin[1][1] );

  pi_div_delta_N[0][0] = PI / ( Nmax[0][0] - Nmin[0][0] );
  pi_div_delta_N[0][1] = PI / ( Nmax[0][1] - Nmin[0][1] );
  pi_div_delta_N[1][0] = pi_div_delta_rc[0][1];
  pi_div_delta_N[1][1] = PI / ( Nmax[1][1] - Nmin[1][1] );

  pi_div_delta_NC[0][0] = PI / ( NCmax[0][0] - Nmin[0][0] );
  pi_div_delta_NC[0][1] = PI / ( NCmax[0][1] - Nmin[0][1] );
  pi_div_delta_NC[1][0] = pi_div_delta_rc[0][1];
  pi_div_delta_NC[1][1] = PI / ( NCmax[1][1] - Nmin[1][1] );

  Q[0][0] = getLineAndConvertToDouble( file );
  Q[0][1] = getLineAndConvertToDouble( file );
  Q[1][0] = Q[0][1];
  Q[1][1] = getLineAndConvertToDouble( file );
  alpha[0][0] = getLineAndConvertToDouble( file );
  alpha[0][1] = getLineAndConvertToDouble( file );
  alpha[1][0] = alpha[0][1];
  alpha[1][1] = getLineAndConvertToDouble( file );
  A[0][0] = getLineAndConvertToDouble( file );
  A[0][1] = getLineAndConvertToDouble( file );
  A[1][0] = A[0][1];
  A[1][1] = getLineAndConvertToDouble( file );

  BIJc_CC1 = getLineAndConvertToDouble( file );
  BIJc_CC2 = getLineAndConvertToDouble( file );
  BIJc_CC3 = getLineAndConvertToDouble( file );
  BIJc_CH1 = getLineAndConvertToDouble( file );
  BIJc_CH2 = getLineAndConvertToDouble( file );
  BIJc_CH3 = getLineAndConvertToDouble( file );
  BIJc_HH1 = getLineAndConvertToDouble( file );
  BIJc_HH2 = getLineAndConvertToDouble( file );
  BIJc_HH3 = getLineAndConvertToDouble( file );

  Beta_CC1 = getLineAndConvertToDouble( file );
  Beta_CC2 = getLineAndConvertToDouble( file );
  Beta_CC3 = getLineAndConvertToDouble( file );
  Beta_CH1 = getLineAndConvertToDouble( file );
  Beta_CH2 = getLineAndConvertToDouble( file );
  Beta_CH3 = getLineAndConvertToDouble( file );
  Beta_HH1 = getLineAndConvertToDouble( file );
  Beta_HH2 = getLineAndConvertToDouble( file );
  Beta_HH3 = getLineAndConvertToDouble( file );

  rho[0][0] = getLineAndConvertToDouble( file );
  rho[0][1] = getLineAndConvertToDouble( file );
  rho[1][0] = rho[0][1];
  rho[1][1] = getLineAndConvertToDouble( file );

  // b) czlon LJ
  rcLJmin[0][0] = getLineAndConvertToDouble( file );
  rcLJmin[0][1] = getLineAndConvertToDouble( file );
  rcLJmin[1][0] = rcLJmin[0][1];
  rcLJmin[1][1] = getLineAndConvertToDouble( file );

  rcLJmax[0][0] = getLineAndConvertToDouble( file );
  rcLJmax[0][1] = getLineAndConvertToDouble( file );
  rcLJmax[1][0] = rcLJmax[0][1];
  rcLJmax[1][1] = getLineAndConvertToDouble( file );

  bLJmin[0][0] = getLineAndConvertToDouble( file );
  bLJmin[0][1] = getLineAndConvertToDouble( file );
  bLJmin[1][0] = bLJmin[0][1];
  bLJmin[1][1] = getLineAndConvertToDouble( file );

  bLJmax[0][0] = getLineAndConvertToDouble( file );
  bLJmax[0][1] = getLineAndConvertToDouble( file );
  bLJmin[1][0] = bLJmin[0][1];
  bLJmax[1][1] = getLineAndConvertToDouble( file );

  rcLJmaxsq[0][0] = rcLJmax[0][0] * rcLJmax[0][0];
  rcLJmaxsq[0][1] = rcLJmax[0][1] * rcLJmax[0][1];
  rcLJmaxsq[1][0] = rcLJmaxsq[0][1];
  rcLJmaxsq[1][1] = rcLJmax[1][1] * rcLJmax[1][1];

  epsilon[0][0] = getLineAndConvertToDouble( file );
  epsilon[0][1] = getLineAndConvertToDouble( file );
  epsilon[1][0] = epsilon[0][1];
  epsilon[1][1] = getLineAndConvertToDouble( file );

  sigma[0][0] = getLineAndConvertToDouble( file );
  sigma[0][1] = getLineAndConvertToDouble( file );
  sigma[1][0] = sigma[0][1];
  sigma[1][1] = getLineAndConvertToDouble( file );

  // c) czlon REBO
  epsilonT_CCCC = getLineAndConvertToDouble( file );
  epsilonT_CCCH = getLineAndConvertToDouble( file );
  epsilonT_HCCH = getLineAndConvertToDouble( file );

  // d) splajny
  // splajny gC1 i gC2
  getLine( file );
  getLine( file );
  getLine( file );
  // liczba wezlow
  number_of_domains = getLineAndConvertToInt( file );
  // gCdom
  for( i = 0; i < number_of_domains; i++ )
    gCdom[i] = getLineAndConvertToDouble( file );
  getLine( file );
  // gC1
  for( i = 0; i < number_of_domains - 1; i++ )
    for( j = 0; j < 6; j++ )
      gC1[i][j] = getLineAndConvertToDouble( file );
  getLine( file );
  // gC2
  for( i = 0; i < number_of_domains - 1; i++ )
    for( j = 0; j < 6; j++ )
      gC2[i][j] = getLineAndConvertToDouble( file );

  // splajn gH
  getLine( file );
  getLine( file );
  getLine( file );
  // liczba wezlow
  number_of_domains = getLineAndConvertToInt( file );
  // gHdom
  for( i = 0; i < number_of_domains; i++ )
    gHdom[i] = getLineAndConvertToDouble( file );
  getLine( file );
  // gH
  for( i = 0; i < number_of_domains - 1; i++ )
    for( j = 0; j < 6; j++ )
      gH[i][j] = getLineAndConvertToDouble( file );

  // splajn pCC
  getLine( file );
  getLine( file );
  getLine( file );
  // liczba wezlow
  number_of_domains = getLineAndConvertToInt( file );
  // pCCdom
  for( i = 0; i < number_of_domains / 2; i++ )
    for( j = 0; j < number_of_domains / 2; j++ )
      pCCdom[i][j] = getLineAndConvertToDouble( file );
  getLine( file );
  // pCC
  for( i = 0; i < (int) pCCdom[0][1]; i++ )
    for( j = 0; j < (int) pCCdom[1][1]; j++ )
      for( k = 0; k < 16; k++ )
        pCC[i][j][k] = getLineAndConvertToDouble( file );

  // splajn pCH
  getLine( file );
  getLine( file );
  getLine( file );
  // liczba wezlow
  number_of_domains = getLineAndConvertToInt( file );
  // pCHdom
  for( i = 0; i < number_of_domains / 2; i++ )
    for( j = 0; j < number_of_domains / 2; j++ )
      pCHdom[i][j] = getLineAndConvertToDouble( file );
  getLine( file );
  // pCH
  for( i = 0; i < (int) pCHdom[0][1]; i++ )
    for( j = 0; j < (int) pCHdom[1][1]; j++ )
      for( k = 0; k < 16; k++ )
        pCH[i][j][k] = getLineAndConvertToDouble( file );

  // splajn piCC
  getLine( file );
  getLine( file );
  getLine( file );
  // liczba wezlow
  number_of_domains = getLineAndConvertToInt( file );
  // piCCdom
  for( i = 0; i < number_of_domains / 2; i++ )
    for( j = 0; j < number_of_domains / 3; j++ )
      piCCdom[i][j] = getLineAndConvertToDouble( file );
  getLine( file );
  // piCC
  for( i = 0; i < (int) piCCdom[0][1]; i++ )
    for( j = 0; j < (int) piCCdom[1][1]; j++ )
      for( k = 0; k < (int) piCCdom[2][1]; k++ )
        for( l = 0; l < 64; l++ )
          piCC[i][j][k][l] = getLineAndConvertToDouble( file );

  // splajn piCH
  getLine( file );
  getLine( file );
  getLine( file );
  number_of_domains = getLineAndConvertToInt( file );
  // piCHdom
  for( i = 0; i < number_of_domains / 2; i++ )
    for( j = 0; j < number_of_domains / 3; j++ )
      piCHdom[i][j] = getLineAndConvertToDouble( file );
  getLine( file );
  // piCH
  for( i = 0; i < (int) piCHdom[0][1]; i++ )
    for( j = 0; j < (int) piCHdom[1][1]; j++ )
      for( k = 0; k < (int) piCHdom[2][1]; k++ )
        for( l = 0; l < 64; l++ )
          piCH[i][j][k][l] = getLineAndConvertToDouble( file );

  // splajn piHH
  getLine( file );
  getLine( file );
  getLine( file );
  number_of_domains = getLineAndConvertToInt( file );
  // piHHdom
  for( i = 0; i < number_of_domains / 2; i++ )
    for( j = 0; j < number_of_domains / 3; j++ )
      piHHdom[i][j] = getLineAndConvertToDouble( file );
  getLine( file );
  // piHH
  for( i = 0; i < (int) piHHdom[0][1]; i++ )
    for( j = 0; j < (int) piHHdom[1][1]; j++ )
      for( k = 0; k < (int) piHHdom[2][1]; k++ )
        for( l = 0; l < 64; l++ )
          piHH[i][j][k][l] = getLineAndConvertToDouble( file );

  // splajn Tij
  getLine( file );
  getLine( file );
  getLine( file );
  number_of_domains = getLineAndConvertToInt( file );
  // Tijdom
  for( i = 0; i < number_of_domains / 2; i++ )
    for( j = 0; j < number_of_domains / 3; j++ )
      Tijdom[i][j] = getLineAndConvertToDouble( file );
  getLine( file );
  // Tijc
  for( i = 0; i < (int) Tijdom[0][1]; i++ )
    for( j = 0; j < (int) Tijdom[1][1]; j++ )
      for( k = 0; k < (int) Tijdom[2][1]; k++ )
        for( l = 0; l < 64; l++ )
          Tijc[i][j][k][l] = getLineAndConvertToDouble( file );

  file.close( );

  thmin = -1.0;
  thmax = -0.995;

  BIJc[0][0][0] = BIJc_CC1;
  BIJc[0][0][1] = BIJc_CC2;
  BIJc[0][0][2] = BIJc_CC3;
  BIJc[0][1][0] = BIJc_CH1;
  BIJc[0][1][1] = BIJc_CH2;
  BIJc[0][1][2] = BIJc_CH3;
  BIJc[1][0][0] = BIJc_CH1;
  BIJc[1][0][1] = BIJc_CH2;
  BIJc[1][0][2] = BIJc_CH3;
  BIJc[1][1][0] = BIJc_HH1;
  BIJc[1][1][1] = BIJc_HH2;
  BIJc[1][1][2] = BIJc_HH3;

  Beta[0][0][0] = Beta_CC1;
  Beta[0][0][1] = Beta_CC2;
  Beta[0][0][2] = Beta_CC3;
  Beta[0][1][0] = Beta_CH1;
  Beta[0][1][1] = Beta_CH2;
  Beta[0][1][2] = Beta_CH3;
  Beta[1][0][0] = Beta_CH1;
  Beta[1][0][1] = Beta_CH2;
  Beta[1][0][2] = Beta_CH3;
  Beta[1][1][0] = Beta_HH1;
  Beta[1][1][1] = Beta_HH2;
  Beta[1][1][2] = Beta_HH3;

  epsilonT[0][0] = epsilonT_CCCC;
  epsilonT[0][1] = epsilonT_CCCH;
  epsilonT[1][0] = epsilonT_CCCH;
  epsilonT[1][1] = epsilonT_HCCH;
}

void AIREBOForceField::allocateMemory( ) {
  cutljsq = new double *[2];
  cutljsq[0] = new double [2];
  cutljsq[1] = new double [2];

  lj1 = new double *[2];
  lj1[0] = new double [2];
  lj1[1] = new double [2];

  lj2 = new double *[2];
  lj2[0] = new double [2];
  lj2[1] = new double [2];

  lj3 = new double *[2];
  lj3[0] = new double [2];
  lj3[1] = new double [2];

  lj4 = new double *[2];
  lj4[0] = new double [2];
  lj4[1] = new double [2];
}

void AIREBOForceField::deallocateMemory( ) {
  delete [] cutljsq[0];
  delete [] cutljsq[1];
  delete [] cutljsq;

  delete [] lj1[0];
  delete [] lj1[1];
  delete [] lj1;

  delete [] lj2[0];
  delete [] lj2[1];
  delete [] lj2;

  delete [] lj3[0];
  delete [] lj3[1];
  delete [] lj3;

  delete [] lj4[0];
  delete [] lj4[1];
  delete [] lj4;

  if ( natoms != 0 )
    deallocateREBO( );
}

void AIREBOForceField::initialize_constants( ) {
  int i, j;
  double tmp_cutmax;
  double tmp;

  for( i = 0; i < 2; i++ )
    for( j = 0; j < 2; j++ ) {
      tmp = pow( sigma[i][j], 6.0 );
      lj4[i][j] = 4.0 * epsilon[i][j] * tmp;
      lj3[i][j] = lj4[i][j] * tmp;
      lj2[i][j] = lj4[i][j] * 6.0;
      lj1[i][j] = lj3[i][j] * 12.0;
      tmp = cutlj * sigma[i][j];
      cutljsq[i][j] = tmp * tmp;
    }

  // promien odciecia dla czlonu REBO zadany przez promien odciecia dla C, gdyz jest on najwiekszy
  cut3rebo = 3.0 * rcMax[0][0];
  cutljrebo = rcLJmax[0][0] + rcMax[0][0];
  cutljrebosq = cutljrebo * cutljrebo;

  cutmax = cut3rebo;

  if ( ljflag != 1 )
    return;

  tmp_cutmax = rcLJmax[0][0] + 2.0 * rcMax[0][0];
  if ( tmp_cutmax > cutmax )
    cutmax = tmp_cutmax;
  tmp_cutmax = cutlj * sigma[0][0];
  if ( tmp_cutmax > cutmax )
    cutmax = tmp_cutmax;
}

void AIREBOForceField::initialize_splines( ) {
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

  // !!!komentarz: make top end of piCC flat instead of zero
  //               also enforces some symmetry
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

  // !!!komentarz: make top end of piCH flat instead of zero
  //               also enforces some symmetry
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
}
