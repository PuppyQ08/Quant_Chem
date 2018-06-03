#include "molec_geo.h"
#include <iostream>
#include <cmath>
#include <armadillo>
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
  _xm = new double[_atomnum];
  _ym = new double[_atomnum];
  _zm = new double[_atomnum];
  _momntInert = new double[3];
  _rotatnCnstnt = new double[3];
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
  delete[] _xm;
  delete[] _ym;
  delete[] _zm;
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

double Molecule::torsionAng(int i, int j, int k, int l){
  double eijk_x = (unit(1,j,i) * unit(2,j,k) - unit(2,j,i) * unit(1,j,k));
  double eijk_y = (unit(2,j,i) * unit(0,j,k) - unit(0,j,i) * unit(2,j,k));
  double eijk_z = (unit(0,j,i) * unit(1,j,k) - unit(1,j,i) * unit(0,j,k));

  double ejkl_x = (unit(1,k,j) * unit(2,k,l) - unit(2,k,j) * unit(1,k,l));
  double ejkl_y = (unit(2,k,j) * unit(0,k,l) - unit(0,k,j) * unit(2,k,l));
  double ejkl_z = (unit(0,k,j) * unit(1,k,l) - unit(1,k,j) * unit(0,k,l));

  double exx = eijk_x * ejkl_x;
  double eyy = eijk_y * ejkl_y;
  double ezz = eijk_z * ejkl_z;

  double tau = (exx + eyy + ezz)/(sin(bondAng(i, j, k))* sin(bondAng(j, k, l)));

  if(tau < -1.0) tau = acos(-1.0);
  else if(tau > 1.0) tau = acos(1.0);
  else tau = acos(tau);

  //to determine the sign of angles
  double cross_x = eijk_y * ejkl_z - eijk_z * ejkl_y;
  double cross_y = eijk_z * ejkl_x - eijk_x * ejkl_z;
  double cross_z = eijk_x * ejkl_y - eijk_y * ejkl_x;
  double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
  cross_x /= norm;
  cross_y /= norm;
  cross_z /= norm;
  double sign = 1.0;
  double dot = cross_x*unit(0,j,k)+cross_y*unit(1,j,k)+cross_z*unit(2,j,k);
  if(dot < 0.0) sign = -1.0;
  return tau*sign;
}

double Molecule::centMass(int judge){
  double sum = 0.0, masum = 0.0, centMass = 0.0;
  if(judge == 0){
    for (int i = 0; i < _atomnum; i++) {
      sum += _x[i] * _amass[_zval[i]];
      masum += _amass[_zval[i]];
    }
    centMass = sum/masum;
    for(int i = 0;i < _atomnum; i++)
      _xm[i] = _x[i]-centMass;//update the coordinate based on cent of mass
    return centMass;
  }
  else if(judge == 1){
    for (int i = 0; i < _atomnum; i++) {
      sum += _y[i] * _amass[_zval[i]];
      masum += _amass[_zval[i]];
    }
    centMass = sum/masum;
    for(int i = 0;i < _atomnum; i++)
      _ym[i] = _y[i]-centMass;
    return centMass;
  }
  else if(judge == 2){
    for (int i = 0; i < _atomnum; i++) {
      sum += _z[i] * _amass[_zval[i]];
      masum += _amass[_zval[i]];
    }
    centMass = sum/masum;
    for(int i = 0;i < _atomnum; i++)
      _zm[i] = _z[i]-centMass;
    return centMass;
  }
  return 0;
}

double* Molecule::momntInert(){
  //to generate inertia tensor matrix
  arma::mat inertia(3,3);
  double Ixx = 0.0, Iyy = 0.0, Izz = 0.0, Ixy = 0.0, Ixz = 0.0, Iyz =0.0;
  for (int i = 0; i < _atomnum; i++) {
    Ixx += _amass[_zval[i]]*(_ym[i]*_ym[i] + _zm[i]*_zm[i]);
    Iyy += _amass[_zval[i]]*(_xm[i]*_xm[i] + _zm[i]*_zm[i]);
    Izz += _amass[_zval[i]]*(_ym[i]*_ym[i] + _xm[i]*_xm[i]);
    Ixy += _amass[_zval[i]]*_xm[i]*_ym[i];
    Ixz += _amass[_zval[i]]*_xm[i]*_zm[i];
    Iyz += _amass[_zval[i]]*_zm[i]*_ym[i];
  }
  inertia(0,0) = Ixx;
  inertia(1,1) = Iyy;
  inertia(2,2) = Izz;
  inertia(0,1) = inertia(1,0) = Ixy;
  inertia(0,2) = inertia(2,0) = Ixz;
  inertia(1,2) = inertia(2,1) = Iyz;
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, inertia);
  _momntInert[0] = eigval(0);
  _momntInert[1] = eigval(1);
  _momntInert[2] = eigval(2);
  return _momntInert;
}

double* Molecule::rotatnCnstnt(){
  double constB = 6.6260755E-34/(8.0*acos(-1.0)*acos(-1.0));
  constB /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
  constB /= 2.99792458E10;;
  _rotatnCnstnt[0] = constB/_momntInert[0];
  _rotatnCnstnt[1] = constB/_momntInert[1];
  _rotatnCnstnt[2] = constB/_momntInert[2];
  return _rotatnCnstnt;
}
