/*
  File:   nlists.h
  Author: vaxquis
 */
#pragma once

#ifndef NLISTS_H
#define	NLISTS_H

namespace AIREBO {

  typedef struct {
    double x, y, z;
    double r, r_sq;
  } vec3d;

  class nLists {
  public:
    int number_of_atoms, *type, *neighbours_num, **neighbours_list;
    vec3d **neighbours_bonds;
  };

}

#endif	/* NLISTS_H */

