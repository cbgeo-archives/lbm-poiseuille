#include "lbgrid.h"

lbgrid::lbgrid(int nx, int ny) : nx_{nx}, ny_{ny} {
  for (unsigned i = 0; i < nx_; ++i) {
    f[i] = new double*[ny_]();
    for (unsigned j = 0; j < ny_; ++j) f[i][j] = new double[Q]();
  }
}

#pragma omp parallel for default(shared) private(i, j) schedule(dynamic, chunk)
void lbgrid::initialize_density(double rho) {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      f[i][j][0] = 4. / 9. * (rho);
      f[i][j][1] = f[i][j][2] = f[i][j][3] = f[i][j][4] = 1. / 9. * (rho);
      f[i][j][5] = f[i][j][6] = f[i][j][7] = f[i][j][8] = 1. / 36. * (rho);
    }
  }
}

double lbgrid::density_function(int i, int j, int k) { return f[i][j][k]; }

#pragma omp parallel for default(shared) private(i) schedule(dynamic, chunk) reducation(+:total_mass)
double lbgrid::sum_density() {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      for (int k = 0; k < Q; ++k) total_mass += f[i][j][k];
    }
  }
  return total_mass;
}

lbgrid::~lbgrid() {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) delete[] f[i][j];
    delete[] f[i];
  }
  delete[] f;
}
