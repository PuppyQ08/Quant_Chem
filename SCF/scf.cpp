#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <string>
#include <unordered_map>
#include "scf.h"
using namespace std;

SCF::SCF(std::string input){
  _nucrepul = 8.002367061810450;
  ifstream iptcoor(input + "coord.txt");
  ifstream iptovlap(input + "ovlap.txt");
  ifstream iptkinetic(input + "kinetic.txt");
  ifstream iptnuclear(input + "nuclear.txt");
  ifstream ipttwoelec(input + "twoelec.txt");
  iptcoor >> _numatom;
  /* normally we need to use Z-value to
  * calculate the total number of orbitals of the molecular
  * but here we just test for water so we basically just type 7 to save time
  */
  _ioff.resize(7);
  _ovlap = new arma::mat(7,7);
  _coreHam = new arma::mat(7,7);
  int i = 0, j = 0,k = 0, l = 0;
  double temp,temp2, temp3;
  /*
  to input overlap integral
  and kinetic
  and nuclear-attraction integral
  */
  while(!iptovlap.eof()){
  iptovlap >> i >> j >> temp;
  iptkinetic >> i >> j >> temp2;
  iptnuclear >> i >> j >> temp3;
  (*_ovlap)(i -1, j -1) = temp;
  (*_ovlap)(j -1, i -1 ) = temp;
  (*_coreHam)(i -1, j -1) = temp2 + temp3;
  (*_coreHam)(j -1, i -1 ) = temp2 + temp3;
}
  /*two electron integral */
  double temp4;
  _ioff[0] = 0;
  for(int i=1; i < 7; i++)
  _ioff[i] = _ioff[i-1] + i;
  while (!ipttwoelec.eof()) {
    ipttwoelec >> i >> j >> k >> l >> temp4;
    _twoelec[getijkl(i, j, k, l)] = temp4;
  }
  iptcoor.close();iptovlap.close();iptkinetic.close();
  ipttwoelec.close();iptnuclear.close();
}

int SCF::getijkl(int i, int j, int k, int l){
  int ij = 0, kl = 0, ijkl = 0;
  ij = (i > j) ? _ioff[i] + j : _ioff[j] + i;
  kl = (k > l) ? _ioff[k] + l : _ioff[k] + l;
  ijkl = (ij > kl) ? _ioff[ij] + kl : _ioff[kl] + ij;
  return ijkl;
}

SCF::~SCF(){
  delete _ovlap;
  delete _coreHam;
}

void SCF::print(){
  for (size_t i = 0; i < 7; i++) {
    printf("%s\n", " ");
    for (size_t j = 0; j < 7; j++) {
    printf("%20.12f, %s", (_OrthogMat)(i, j), " ");
    }
  }
}

void SCF::calculation(){
  /*diagonlize overlap Matrix*/
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, *_ovlap);
  //to get diagonlized eigenvalue matrix
  //arma::mat eigvalmat = eigvec.t() * (*_ovlap) * eigvec;
  // to get Orthogonalization Matrix
  //_OrthogMat = eigvec * (arma::sqrtmat_sympd(eigvalmat).i()) * eigvec.t();

}
