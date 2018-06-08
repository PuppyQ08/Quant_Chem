#ifndef HARMNC_H
#define HARMNC_H
#include <iostream>
#include <fstream>
#include <armadillo>
class Harmnc_vib{
public:
  /*during the constructor
  *we would have input zvalue, coordinates and hessian function
  */
  Harmnc_vib();
  ~Harmnc_vib();
  void printout();
  /*to diagonlize mass-weighted hessian matrix*/
  void freq();
private:
  size_t _numatom;
  double **_coor;
  double **_hessian;
  int *_zval;
  double _amass[9]= {0.0, 1.00782503223,  4.00260325413,  7.0160034366,  9.012183065, 11.00930536, 12.0000000, 14.00307400443, 15.99491461957};
  arma::mat *_msHessian;
};

#endif
