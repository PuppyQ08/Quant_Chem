#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <string>
#include "scf.h"
using namespace std;

SCF::SCF(std::string input){
  ifstream iptcoor(input + "coord.txt");
  ifstream iptovlap(input + "ovlap.txt");
  iptcoor >> _numatom;
  /* normally we need to use Z-value to calculate the total number of orbitals of the molecular
  * but here we just test for water so we basically just type 7 to save time
  */
  _ovlap = new arma::mat(7,7);
  int i = 0, j = 0;
  double temp;
  cout<<"here"<<endl;
  while(!iptovlap.eof()){
  iptovlap >> i >> j >> temp;
  //cout<< i<<endl;
  (*_ovlap)(i -1, j -1) = temp;
  (*_ovlap)(j -1, i -1 ) = temp;
}
}

SCF::~SCF(){
  delete _ovlap;
}

void SCF::print(){
  for (size_t i = 0; i < 7; i++) {
    printf("%s\n", " ");
    for (size_t j = 0; j < 7; j++) {
    printf("%20.12f, %s", (*_ovlap)(i, j), " ");
    }
  }
}
