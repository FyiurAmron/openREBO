/*
    File:   main.cpp
    Author: vaxquis
 */

#include <cstdlib>
#include <cstdio>
#include <ctime>

#include "OpenREBO.h"

using namespace std;
using namespace OpenREBO;

int tmp, total_read, total_compute;

void OpenREBO_test( ) {
  clock_t t1;
  printf( "AIREBO test...\n" );
  AIREBO aff( "CH.airebo", 0.0, false, false, 256 );
  //aff.readNList( "small.nlists" );
  t1 = clock( );
  aff.readNList( "big.nlists" );
  tmp = clock( ) - t1;
  cout << "read() clocks used: " << tmp << endl;
  total_read += tmp;
  t1 = clock( );
  aff.compute( );
  tmp = clock( ) - t1;
  cout << "compute() clocks used: " << tmp << endl;
  total_compute += tmp;
  cout << aff.getEnergyTotal( ) << endl;
}

int main( int argc, char** argv ) {
  setbuf( stdout, NULL );
  ios::sync_with_stdio( false );
  for( int i = 0; i < 10; i++ )
    OpenREBO_test( );
  cout << "total_read: " << total_read << endl
          << "total_compute: " << total_compute << endl;

  return 0;
}

