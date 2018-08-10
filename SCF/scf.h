#ifndef SCF_H
#define SCF_H
#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>
#include <unordered_map>
class SCF{
public:
  //input constructor for reading file
  SCF(std::string input);
  ~SCF();//destructor
  void calculation();//doing scf calculation
  void print();//fn to print out the result
  //helper fn to get ijkl for two electron integral
  int getijkl(int i, int j, int k, int l);
private:
  int _numatom;
  double _nucrepul;
  arma::mat * _ovlap;
  arma::mat * _coreHam;
  std::unordered_map<int, double> _twoelec;
  std::vector<int> _ioff;//this one is for storing i(i+1)/2
  arma::mat _OrthogMat;
};

#endif
