/*
    File:   vector.h
    Author: vaxquis
 */
#pragma once

#ifndef VECTOR_H
#define	VECTOR_H

#include <cmath>

template <typename T>
void sum( T v0[3], const T v1[3], const T v2[3] ) {
  v0[0] = v1[0] + v2[0];
  v0[1] = v1[1] + v2[1];
  v0[2] = v1[2] + v2[2];
}

template <typename T>
void diff( T v0[3], const T v1[3], const T v2[3] ) {
  v0[0] = v1[0] - v2[0];
  v0[1] = v1[1] - v2[1];
  v0[2] = v1[2] - v2[2];
}

template <typename T>
void neg( T v0[3], const T v1[3] ) {
  v0[0] = -v1[0];
  v0[1] = -v1[1];
  v0[2] = -v1[2];
}

template <typename T>
void scale( T v0[3], const T v1[3], const T scale ) {
  v0[0] = scale * v1[0];
  v0[1] = scale * v1[1];
  v0[2] = scale * v1[2];
}

template <typename T>
T dot( const T v0[3], const T v1[3] ) {
  return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
}

template <typename T>
void cross( T v0[3], const T v1[3], const T v2[3] ) {
  v0[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v0[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v0[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

template <typename T >
T length_sq( const T v[3] ) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

template <typename T >
T length( const T v[3] ) {
  return sqrt( length_sq( v ) );
}

template <typename T >
T cos_theta_clamp( const T v0[3], const T v1[3], T len0, T len1 ) {
  T val = ( v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] ) / ( len0 * len1 );
  if ( val > 1.0 )
    return 1.0;
  else if ( val < 1.0 )
    return 1.0;
  return val;
}

template <typename T >
T cos_theta_clamp( const T v0[3], const T v1[3], T scale ) {
  T val = ( v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] ) * scale;
  if ( val > 1.0 )
    return 1.0;
  else if ( val < 1.0 )
    return 1.0;
  return val;
}

#endif	/* VECTOR_H */

