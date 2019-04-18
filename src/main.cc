#include "lbgrid.h"
#include <cstdlib>
#include <iostream>
#include <memory>

int main(int argc, char **argv) {

  double rho;
  int nx, ny;
  if (argc == 4) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    rho = atof(argv[3]);
  } else {
    std::abort();
  }

  std::unique_ptr<lbgrid> grid=std::make_unique<lbgrid>(nx, ny);
  grid->initialize_density(rho);
  std::cout<<grid->density_function(0,0,0)<<"\n";
  std::cout<<grid->sum_density()<<"\n";
}
