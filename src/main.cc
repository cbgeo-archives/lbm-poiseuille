#include "lbgrid.h"
#include <cstdlib>
#include <memory>

int main(int argc, char** argv) {

  double rho_0, Fx, Fy;
  int nx, ny;
  if (argc == 6) {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    rho_0 = atof(argv[3]);
    Fx = atof(argv[4]);
    Fy = atof(argv[5]);
  } else {
    std::abort();
  }

  double tau = 0.933;
  double dt = 1.;

  std::unique_ptr<lbgrid> grid = std::make_unique<lbgrid>(nx, ny);
  grid->initialize_density(rho_0);
  std::cout << grid->f_(0, 0, 0) << "\n";
  std::cout << grid->sum_density() << "\n";

  for (int i = 0; i < 3; ++i) {
    std::cout << i << std::endl;
    grid->compute_macro_var(dt,Fx,Fy);
    std::cout << grid->ux_(0, 0) << "\n";
    std::cout << grid->uy_(0, 0) << "\n";
    grid->equilibrium_density(dt);
    std::cout << grid->feq_(0, 0, 0) << "\n";
    grid->collision(tau, dt, Fx, Fy);
    std::cout << grid->fcol_(0, 0, 0) << "\n";
    grid->streaming();
    std::cout << grid->f_(0, 0, 0) << "\n";
    std::cout << grid->sum_density() << "\n";
  }

  std::cout<<"Ux at x=0"<<std::endl;
for (int j=20;j<ny;++j) {
  std::cout<< grid->uxsection(0)[j]<<std::endl;
}
}
