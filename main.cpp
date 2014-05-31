/*
 * File:   main.cpp
 * Author: toor
 *
 * Created on 31 maja 2014, 20:51
 */

#include <cstdlib>
#include <cstdio>

#include "airebo_force_field.h"

using namespace std;

/*
 *
 */
int main( int argc, char** argv ) {
  printf( "test\n" );
  AIREBOForceField* aff = new AIREBOForceField( "CH.airebo", 1.0, false, false, 256 );
  //aff->compute()
  return 0;
}

