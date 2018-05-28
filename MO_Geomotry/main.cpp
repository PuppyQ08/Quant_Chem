#include <iostream>
#include <fstream>
//#include <iomanip>
#include <string>
#include "molec_geo.h"
 int main(int argc, char const *argv[]) {
  if(argc == 2){
    if(std::string(argv[1]) == "acetal"){
      Molecule obj("acetaldehyde.txt", 0);
      obj.printfn();
    }
  }

  return 0;
}
