#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <string>
#include "scf.h"
using namespace std;

SCF::SCF(std::string input){
  _nucrepul = 8.002367061810450;
  ifstream iptcoor(input + "coord.txt");
  ifstream iptovlap(input + "ovlap.txt");
  ifstream iptkinetic(input + "kinetic.txt");
  ifstream iptnuclear(input + "nuclear.txt");
  iptcoor >> _numatom;
  /* normally we need to use Z-value to calculate the total number of orbitals of the molecular
  * but here we just test for water so we basically just type 7 to save time
  */
  _ovlap = new arma::mat(7,7);
  _coreHam = new arma::mat(7,7);
  int i = 0, j = 0;
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
}

SCF::~SCF(){
  delete _ovlap;
  delete _coreHam;
}

void SCF::print(){
  for (size_t i = 0; i < 7; i++) {
    printf("%s\n", " ");
    for (size_t j = 0; j < 7; j++) {
    printf("%20.12f, %s", (*_coreHam)(i, j), " ");
    }
  }
}
