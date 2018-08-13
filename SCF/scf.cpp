#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <string>
#include <unordered_map>
#include "scf.h"
using namespace std;
using namespace arma;
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
  * but here we just test for water so we basically just type _numorbit to save time
  */
  _numorbit = 7;
  _numoccp = 5;
  _ioff.resize(pow(_numorbit,4));
  //to set up ioff:
  _ioff[0] = 0;
  for(int i=1; i < pow(_numorbit,4); i++)
  _ioff[i] = _ioff[i-1] + i;
  //
  _ovlap = new arma::mat(_numorbit,_numorbit);
  _coreHam = new arma::mat(_numorbit,_numorbit);
  _Pini = new arma::mat(_numorbit,_numorbit);
  _Pnext = new arma::mat(_numorbit,_numorbit);
  _Fini = new arma::mat(_numorbit,_numorbit);
  _Fnext = new arma::mat(_numorbit,_numorbit);
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
  kl = (k > l) ? _ioff[k] + l : _ioff[l] + k;
  ijkl = (ij > kl) ? _ioff[ij] + kl : _ioff[kl] + ij;
  return ijkl;
}

SCF::~SCF(){
  delete _ovlap;
  delete _coreHam;
}

void SCF::print(arma::mat ipt){
  for (size_t i = 0; i < _numorbit; i++) {
    printf("%s\n", " ");
    for (size_t j = 0; j < _numorbit; j++) {
    printf("%20.12f, %s", (ipt)(i, j), " ");
    }
  }
}


void SCF::calculation(){
  /*diagonlize overlap Matrix*/
  arma::vec Seigval;
  arma::mat Seigvec;
  arma::eig_sym(Seigval, Seigvec, *_ovlap);
  /*to get diagonlized eigenvalue matrix
  *arma::mat eigvalmat = eigvec.t() * (*_ovlap) * eigvec;//works good
  but the following one would be easier to use*/
  arma::mat Seigvalmat = arma::diagmat(Seigval);
  // to get Orthogonalization Matrix
  //So element-wise inverse and square-root ! except ij term in matrix!

  /*arma::mat et = 1/ sqrt(abs(eigvalmat));
  *arma::mat S_ihalf = arma::eye(_numorbit,_numorbit);
  for (size_t i = 0; i < _numorbit; i++) {
      S_ihalf(i, i) = sqrt(1.0 /eigvalmat(i,i));
  } this one works but the following one would be more neat
  */
  arma::mat lmd_sqrti = arma::sqrt(arma::inv(Seigvalmat));

  _Ssqrtinv = Seigvec * lmd_sqrti * Seigvec.t();//different from website But I believe it caused by different order of eigenvalue? No it is caused by wrong input file
  /*temporaly moving on
  * get the Fock matrix for inital guess from hailtonian
  * then transfer to orthogonal basis to get Fock prime
  * then diagonlize Fock! in the orthogonal basis!
  * then transfer back to original basis which is the atomic orbitals
  */
  *_Fini = *_coreHam;
  _Finiprime = _Ssqrtinv.t() * (*_coreHam) * _Ssqrtinv;
  //digonalize Fock Matrix
  arma::vec Feigval;
  arma::mat Feigvec;
  arma::eig_sym(Feigval, Feigvec, _Finiprime);
  arma::mat Coe_orig = _Ssqrtinv * Feigvec;
  //density matrix equals to C * C.t()
  *_Pini = Coe_orig.cols(0, _numoccp -1) *  Coe_orig.cols(0, _numoccp -1).t();
  //then get the Energy
  arma::mat Eini = (*_Pini) % ((*_coreHam) + *_Fini);
  _Eini = arma::accu(Eini);
  _Eini += _nucrepul;
  printf("%f\n", _Eini);

  /*iteration!
  */
  while(true){
  /*build new Fock ! by old density matrix */
  for (size_t i = 0; i < _numorbit; i++) {
      for (size_t j = 0; j < _numorbit; j++) {
        (*_Fnext)(i,j) = (*_coreHam)(i,j);
        for (size_t k = 0; k < _numorbit; k++) {
          for (size_t l = 0; l < _numorbit; l++) {
            (*_Fnext)(i,j) += (*_Pini)(k,l) * (2 * _twoelec[getijkl(i+1,j+1,k+1,l+1)] - _twoelec[getijkl(i+1,k+1,j+1,l+1)]);
          }
        }
      }
    }
    _Finiprime = _Ssqrtinv.t() * (*_Fnext) * _Ssqrtinv;
    /* to get density matrix and get Energy again..
    */
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, _Finiprime);
    arma::mat Coe = _Ssqrtinv * eigvec;
    //density matrix equals to C * C.t()
    *_Pnext = Coe.cols(0, _numoccp -1) *  Coe.cols(0, _numoccp -1).t();
    //print(*_Pnext);
    _Enext = arma::accu((*_Pnext) % ((*_coreHam) + *_Fnext)) + _nucrepul;
    printf("%20.12f\n", _Enext);
    if(abs(_Enext - _Eini) < 0.00000000001)
    break;
    _Eini = _Enext;
    *_Fini = *_Fnext;
    *_Pini = *_Pnext;
  }
}
