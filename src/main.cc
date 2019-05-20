#include "lbgrid.h"
#include <cstdlib>
#include <iomanip>
#include <memory>

// #define _log_

int main(int argc, char** argv) {

  double rho_0, Fx, Fy;
  int nx, ny;
  double nstep = 0.;
  if (argc == 7) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    rho_0 = atof(argv[3]);
    Fx = atof(argv[4]);
    Fy = atof(argv[5]);
    nstep = atof(argv[6]);
  } else {
    std::abort();
  }

  double tau = 0.933;

  std::unique_ptr<lbgrid> grid = std::make_unique<lbgrid>(nx, ny);
  grid->initialize_density(rho_0);
  std::cout << grid->f_(2, 2, 1) << "\n";
  std::cout << grid->sum_density() << "\n";

  for (int nt = 0; nt <= nstep; ++nt) {
#ifdef _log_
    std::cout << "step number = " << nt << std::endl;
#endif
    grid->compute_macro_var(Fx, Fy);
    grid->equilibrium_density();
    grid->collision(tau, Fx, Fy);
    grid->streaming();
    grid->analytical_solution(tau, nt, Fx);
    grid->write_vtk(nt);
    grid->compute_macro_var(Fx, Fy);
  }

  std::cout << "ux at x=2" << std::endl;
  for (int j = 0; j < ny; ++j) {
    std::cout << grid->uxsection(2)[j] << std::endl;
  }
  std::cout << "u_an at x=2" << std::endl;
  for (int j = 0; j < ny; ++j) {
    std::cout << std::setprecision(10) << grid->u_an_section(2)[j] << std::endl;
  }
}
