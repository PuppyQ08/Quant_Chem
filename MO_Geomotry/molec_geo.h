#ifndef MOLEC_GEO_H
#define MOLEC_GEO_H

#include <iostream>
#include <fstream>
//#include <iomanip>

class Molecule{
public:
  /*constructor and destructor
  *since I wont use overloaded= in any case,
  *I may leave that one out from rule of three
  */
  Molecule(const char *filename, int charge);
  ~Molecule();
  /*varies functions*/
  /*printfn*/
  void printfn();

  /*bond distance calculation*/
  void bondLength();

private:
  int _atomnum;
  int _charge;
  int *_zval;
  double *_x;
  double *_y;
  double *_z;
  double **_distc;
};
#endif
