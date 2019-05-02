#include "lbgrid.h"
#include <cstdlib>
#include <memory>

int main(int argc, char** argv) {

  double rho_0, Fx, Fy;
  int nx, ny, nstep;
  if (argc == 7) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    rho_0 = atof(argv[3]);
    Fx = atof(argv[4]);
    Fy = atof(argv[5]);
    nstep = atoi(argv[6]);
  } else {
    std::abort();
  }

  double tau = 0.933;
  double dt = 1.;

  std::unique_ptr<lbgrid> grid = std::make_unique<lbgrid>(nx, ny);
  grid->initialize_density(rho_0);
  std::cout << grid->f_(2, 2, 1) << "\n";
  std::cout << grid->sum_density() << "\n";

  for (int i = 0; i <= nstep; ++i) {
     std::cout << "step number = "<<i << std::endl;
    grid->compute_macro_var(dt, Fx, Fy);
    std::cout << "ux(2,2) = " << grid->ux_(2, 2) << "\n";
    std::cout << "uy(2,2) = " << grid->uy_(2, 2) << "\n";
    grid->equilibrium_density(dt);
     std::cout << "feq(2,2,1) = "<<grid->feq_(2, 2, 1) << "\n";
    grid->collision(tau, dt, Fx, Fy);
     std::cout << "fcol(2,2,1) = "<<grid->fcol_(2, 2, 1) << "\n";
    grid->streaming();
     std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 1) << "\n";
     std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 2) << "\n";
     std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 3) << "\n";
     std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 4) << "\n";
     std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 5) << "\n";
     std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 6) << "\n";
     std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 7) << "\n";
     std::cout << "f(2,2,1) = "<<grid->f_(2, 2, 8) << "\n";
     std::cout << "sum_density = "<<grid->sum_density() << "\n";
    grid->compute_macro_var(dt,Fx,Fy);
    grid->write_vtk(i);
  }

  std::cout << "Ux at x=2" << std::endl;
  for (int j = 0; j < ny; ++j) {
    std::cout << grid->uxsection(2)[j] << std::endl;
  }
}
