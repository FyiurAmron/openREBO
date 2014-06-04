/*
    File:   main.cpp
    Author: vaxquis
 */

#include <cstdlib>
#include <cstdio>

#include "OpenREBO.h"

using namespace std;
using namespace OpenREBO;

int main( int argc, char** argv ) {
  printf( "test\n" );
  AIREBO aff( "CH.airebo", 0.0, false, false, 256 );
  aff.readNList( "small.nlists" );
  aff.compute( );
  cout << aff.getEnergyTotal( ) << endl;
  return 0;
}

