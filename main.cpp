/*
 * File:   main.cpp
 * Author: toor
 *
 * Created on 31 maja 2014, 20:51
 */

#include <cstdlib>
#include <cstdio>

#include "airebo_force_field.h"
#include "nlist.h"

using namespace std;
using namespace AIREBO;

/*
 *
 */
int main( int argc, char** argv ) {
  printf( "test\n" );
  ForceField aff( "CH.airebo", 0.0, false, false, 256 );
  aff.readNList( "small.nlists" );
  aff.compute( );
  cout << aff.getEnergyTotal();
  return 0;
}

