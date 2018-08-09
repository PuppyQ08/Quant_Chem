#ifndef SCF_H
#define SCF_H
#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>
class SCF{
public:
  //input constructor for reading file
  SCF(std::string input);
  ~SCF();//destructor
  //void calculation();//doing scf calculation
  void print();//fn to print out the result
private:
  int _numatom;
  double _nucrepul;
  arma::mat * _ovlap;
};

#endif
