#include "lbgrid.h"
#include <cstdlib>
#include <iostream>

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

  lbgrid grid(nx, ny);
  grid.initialize_density(rho);
  grid.print_f();
}
