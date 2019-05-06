#include "lbgrid.h"
#include <cstdlib>
#include <iomanip>
#include <memory>

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

  double tau = 0.51;
  double dt = 1.;

  std::unique_ptr<lbgrid> grid = std::make_unique<lbgrid>(nx, ny);
  grid->initialize_density(rho_0);
  std::cout << grid->f_(2, 2, 1) << "\n";
  std::cout << grid->sum_density() << "\n";

  for (int nt = 0; nt <= nstep; ++nt) {
    std::cout << "step number = " << nt << std::endl;
    grid->compute_macro_var(dt, Fx, Fy);
    // std::cout << "ux(2,2) = " << grid->ux_(2, 2) << "\n";
    // std::cout << "uy(2,2) = " << grid->uy_(2, 2) << "\n";
    grid->equilibrium_density(dt);
    //  std::cout << "feq(2,2,1) = "<<grid->feq_(2, 2, 1) << "\n";
    grid->collision(tau, dt, Fx, Fy);
    //  std::cout << "fcol(2,2,1) = "<<grid->fcol_(2, 2, 1) << "\n";
    grid->streaming();
    /* std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 1) << "\n";
     std::cout << "f(2,2,2) = "<<grid->f_(2, 2, 2) << "\n";
     std::cout << "f(2,2,3) = "<<grid->f_(2, 2, 3) << "\n";
     std::cout << "f(2,2,4) = "<<grid->f_(2, 2, 4) << "\n";
     std::cout << "f(2,2,5) = "<<grid->f_(2, 2, 5) << "\n";
     std::cout << "f(2,2,6) = "<<grid->f_(2, 2, 6) << "\n";
     std::cout << "f(2,2,7) = "<<grid->f_(2, 2, 7) << "\n";
     std::cout << "f(2,2,8) = "<<grid->f_(2, 2, 8) << "\n";
     std::cout << "sum_density = "<<grid->sum_density() << "\n";*/
    grid->compute_macro_var(dt, Fx, Fy);
    // grid->analytical_solution(tau, nt);
    if (nt % 10000 == 0) grid->write_vtk(nt);
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
