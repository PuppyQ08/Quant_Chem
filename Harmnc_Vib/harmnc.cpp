#include <iostream>
#include <fstream>
#include <armadillo>
#include "harmnc.h"

using namespace std;
Harmnc_vib::Harmnc_vib(){
  //create dynamic memory
  ifstream iptcoor("benzene.txt");
  ifstream iptHessian("benzHessian.txt");
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
      (*_msHessian)(3*i,3*j) = _hessian[3*i][3*j]/(_amass[_zval[i]]*_amass[_zval[j]]);
      (*_msHessian)(3*i + 1,3*j + 1) = _hessian[3*i + 1][3*j + 1]/(_amass[_zval[i]]*_amass[_zval[j]]);
      (*_msHessian)(3*i + 2,3*j + 2) = _hessian[3*i + 2][3*j + 2]/(_amass[_zval[i]]*_amass[_zval[j]]);
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
}

void Harmnc_vib::printout(){
  for (size_t i = 0; i < 3*_numatom; i++) {
    //for (size_t j = 0; j < _numatom; j++) {
    printf("%20.12f\n", (*_msHessian)(0,i));
  //}
  }

}
