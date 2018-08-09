#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include "harmnc.h"

using namespace std;
Harmnc_vib::Harmnc_vib(){
  //create dynamic memory
  ifstream iptcoor("water.txt");
  ifstream iptHessian("waterHessian.txt");
  iptcoor >> _numatom;
  iptHessian >> _numatom;
  //coordinates memory
  _coor = new double*[_numatom];
  for (size_t i = 0; i < _numatom; i++) {
    _coor[i] = new double[3];
  }
  _zval = new int[_numatom];
  //hessian memory
  _hessian = new double*[3*_numatom];
  for (size_t i = 0; i < 3*_numatom; i++) {
    _hessian[i] = new double[3*_numatom];
  }
  //read the file data
  for (size_t i = 0;i < _numatom; i++) {
    iptcoor>> _zval[i] >>_coor[i][0] >> _coor[i][1] >> _coor[i][2];
  }
  for (size_t i = 0; i < 3*_numatom; i++) {
    for (size_t j = 0; j < 3*_numatom; j++) {
      iptHessian >> _hessian[i][j];
    }
  }
  //to create mass-weighted Hessian Matrix
  _msHessian = new arma::mat(3*_numatom, 3*_numatom);
  for (size_t i = 0; i < _numatom; i++) {
    for (size_t j = 0; j < _numatom; j++) {
      for (size_t k = 0; k < 3; k++) {
        for (size_t l = 0; l < 3; l++) {
      (*_msHessian)(3*i + k,3*j + l) = _hessian[3*i + k][3*j + l]/sqrt(_amass[_zval[i]]*_amass[_zval[j]]);
    }
    }
  }
}
}

Harmnc_vib::~Harmnc_vib(){
  for (size_t i = 0; i < _numatom; i++) {
    delete _coor[i];
  }
  delete []_coor;
  for (size_t i = 0; i < 3*_numatom; i++) {
    delete _hessian[i];
  }
  delete []_hessian;
  delete _msHessian;// I am not sure whether if this work or not
}

void Harmnc_vib::printout(){
  for (size_t i = 0; i < 3*_numatom; i++) {
    //for (size_t j = 0; j < _numatom; j++) {
    printf("%20.12f\n", (*_msHessian)(0,i));
  //}
  }
}
void Harmnc_vib::freq(){
  arma::mat eigvec;
  arma::vec eigval;
  arma::eig_sym(eigval, eigvec, (*_msHessian));
  double constv =4.359743E-18/(1.660538E-27 * 0.5292E-10 * 0.5292E-10);
  double *freq =new double[3*_numatom];
  for (size_t i = 0; i < 3*_numatom; i++) {
    //printf("%20.12f\n",eigval(i));
    freq[i] =sqrt(abs(eigval(i)) * constv)/(2.99792458E10 * 2 * acos(-1.0)); //remember to include 2 pi here!! w = 2pi* v
    printf("%20.12f\n", freq[i]);
  }
  delete []freq;
}
