#include "molec_geo.h"
#include <iostream>
#include <cmath>
using namespace std;
Molecule::Molecule(const char *filename, int charge){
  _charge = charge;
  ifstream iptfile(filename);
  iptfile >> _atomnum;
  natom = _atomnum;
  _zval = new int[_atomnum];
  _x = new double[_atomnum];
  _y = new double[_atomnum];
  _z = new double[_atomnum];
  //  _distc = new double*[_atomnum];
  for(int i = 0; i < _atomnum; i++){
    iptfile>> _zval[i] >> _x[i]>> _y[i]>> _z[i];
    //  _distc[i] = new double[_atomnum];
    }
}
Molecule::~Molecule(){
  delete[] _zval;
  delete[] _x;
  delete[] _y;
  delete[] _z;
}

double Molecule::unit(int judge, int i, int j){
  if(judge == 0)
    return -(_x[i] - _x[j])/bondLength(i,j);
    if(judge == 1)
      return -(_y[i] - _y[j])/bondLength(i,j);
      if(judge == 2)
        return -(_z[i] - _z[j])/bondLength(i,j);
    return 0;
}

void Molecule::printfn(){
  for(int i = 0; i < _atomnum; i++){
    printf(" %d %20.12f %20.12f %20.12f\n", _zval[i], _x[i], _y[i], _z[i]);
  }
}

double Molecule::bondLength(int i, int j){
  double bondL= sqrt(pow((_x[i] - _x[j]), 2) + pow((_y[i]- _y[j]), 2) + pow((_z[i] - _z[j]),2));
      //printf("%20.12f %d-%d\n",bondL, i ,j);
  return bondL;
}

double Molecule::bondAng(int i, int j, int k){

    //double bondA = acos( ((_x[j] - _x[i])*(_x[j] - _x[k]) + (_y[j] - _y[i])*(_y[j] - _y[k]) + (_z[j] - _z[i])* (_z[j] - _z[k]))/(bondLength(i,j) * bondLength(j,k)));
  //  double bondA = acos(bondphase);
    double bondA = acos(unit(0,j,i) * unit(0,j,k) + unit(1,j,i) * unit(1,j,k) + unit(2,j,i) * unit(2,j,k));
    return bondA;
}

double Molecule::bondplanA(int i, int j,int k, int l){
    double ejkl_x = (unit(1,k,j) * unit(2,k,l) - unit(2,k,j) * unit(1,k,l));
    double ejkl_y = (unit(2,k,j) * unit(0,k,l) - unit(0,k,j) * unit(2,k,l));
    double ejkl_z = (unit(0,k,j) * unit(1,k,l) - unit(1,k,j) * unit(0,k,l));

    double exx = ejkl_x * unit(0,k,i);
    double eyy = ejkl_y * unit(1,k,i);
    double ezz = ejkl_z * unit(2,k,i);
    double theta = (exx + eyy + ezz)/sin(bondAng(j,k,l));

    if(theta < -1.0) theta = asin(-1.0);
    else if(theta > 1.0) theta = asin(1.0);
    else theta = asin(theta);

    return theta;
}

double torsionAng(int i, int j, int k, int l){
  
}
