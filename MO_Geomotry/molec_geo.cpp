#include "molec_geo.h"
#include <iostream>
using namespace std;
Molecule::Molecule(const char *filename, int charge){
    _charge = charge;
    ifstream iptfile(filename);
    iptfile >> _atomnum;
    _zval = new int[_atomnum];
    _x = new double[_atomnum];
    _y = new double[_atomnum];
    _z = new double[_atomnum];
    _distc = new double*[_atomnum];
    for(int i = 0; i < _atomnum; i++){
      iptfile>> _zval[i] >> _x[i]>> _y[i]>> _z[i];
      _distc[i] = new double[_atomnum];
    }
}
Molecule::~Molecule(){
      delete[] _zval;
      delete[] _x;
      delete[] _y;
      delete[] _z;
}

void Molecule::printfn(){
  for(int i = 0; i < _atomnum; i++){
       printf("%d %20.12f %20.12f %20.12f\n", (int) _zval[i], _x[i], _y[i], _z[i]);
  }
}

void Molecule::bondLength(){
  for(int i = 0; i < _atomnum; i++){
    for(int j = 0; j < _atomnum; j++){
      
    }
  }
}
