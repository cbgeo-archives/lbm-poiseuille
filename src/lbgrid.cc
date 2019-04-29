#include "lbgrid.h"

lbgrid::lbgrid(int nx, int ny) : nx_{nx}, ny_{ny} {
  for (unsigned i = 0; i < nx_; ++i) {
    f[i] = new double*[ny_]();
    ux[i] = new double[ny_]();
    uy[i] = new double[ny_]();
    feq[i] = new double*[ny_]();
    fcol[i] = new double*[ny_]();
    for (unsigned j = 0; j < ny_; ++j) {
      f[i][j] = new double[Q]();
      feq[i][j]=new double[Q]();
      fcol[i][j]=new double[Q]();
    }
  }
}

#pragma omp parallel for default(shared) private(i, j) schedule(dynamic, chunk)
void lbgrid::initialize_density(double rho_0) {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      f[i][j][0] = w0 * (rho_0);
      f[i][j][1] = f[i][j][2] = f[i][j][3] = f[i][j][4] = w1_4 * (rho_0);
      f[i][j][5] = f[i][j][6] = f[i][j][7] = f[i][j][8] = w5_8 * (rho_0);
    }
  }
}

#pragma omp parallel for default(shared) private(i) schedule(dynamic, chunk) reducation(+:total_mass)
double lbgrid::sum_density() {
  double total_mass = 0.0;
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      for (int k = 0; k < Q; ++k) total_mass += f[i][j][k];
    }
  }
  return total_mass;
}

void lbgrid::compute_macro_var() {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      rho = f[i][j][0]+f[i][j][1] + f[i][j][2] + f[i][j][3] + f[i][j][4] +
                  f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];
      ux[i][j] = (f[i][j][1] + f[i][j][5] + f[i][j][8] -
                  (f[i][j][3] + f[i][j][6] + f[i][j][7])) /
                 rho;
      uy[i][j] = (f[i][j][2] + f[i][j][5] + f[i][j][6] -
                  (f[i][j][4] + f[i][j][7] + f[i][j][8])) /
                 rho;
    }
  }
}

void lbgrid::equilibrium_density(double Fx, double Fy, double rho_0, const double dt) {

  double a1 = 3.;
  double a2 = 9. / 2.;
  double a3 = 3. / 2.;

  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {

      const double uxeq = ux[i][j] + dt * Fx / (2*rho_0);
      const double uyeq = uy[i][j] + dt * Fy / (2*rho_0);
      const double uxeqsq = uxeq* uxeq;
      const double uyeqsq = uyeq * uyeq;
      const double ueqsq = uxeqsq + uyeqsq;

      rho = f[i][j][0]+f[i][j][1] + f[i][j][2] + f[i][j][3] + f[i][j][4] +
                  f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];

      feq[i][j][0] = w0 * rho * (1. - a3 * ueqsq);
      feq[i][j][1] =
          w1_4 * rho *
          (1. + a1 * uxeq + a2 * uxeqsq - a3 * ueqsq);
      feq[i][j][2] =
          w1_4 * rho *
          (1. + a1 * uyeq + a2 * uyeqsq - a3 * ueqsq);
      feq[i][j][3] =
          w1_4 * rho*
          (1. - a1 * uxeq + a2 * uxeqsq - a3 * ueqsq);
      feq[i][j][4] =
          w1_4 * rho *
          (1. - a1 * uyeq + a2 * uyeqsq - a3 * ueqsq);
      feq[i][j][5] = w5_8 * rho *
                     (1. + a1 * (uxeq + uyeq) +
                      a2 * (uxeqsq + uyeqsq) - a3 * ueqsq);
      feq[i][j][6] = w5_8 * rho *
                     (1. + a1 * (-uxeq + uyeq) +
                      a2 * (-uxeqsq + uyeqsq) - a3 * ueqsq);
      feq[i][j][7] = w5_8 * rho *
                     (1. + a1 * (-uxeq - uyeq) +
                      a2 * (-uxeqsq - uyeqsq) - a3 * ueqsq);
      feq[i][j][8] = w5_8 * rho *
                     (1. + a1 * (uxeq - uyeq) +
                      a2 * (uxeqsq - uyeqsq) - a3 * ueqsq);
    }
  }
}

void lbgrid::collision(double tau) {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      for (int q = 0; q < Q; ++q) {
        fcol[i][j][q] = f[i][j][q] - (f[i][j][q] - feq[i][j][q]) / tau;
      }
    }
  }
}

void lbgrid::streaming() {
  for (int i = 0; i < nx_; ++i) {
    int in = (i > 0) ? (i - 1) : (nx_ - 1);
    int ip = (i < nx_ - 1) ? (i + 1) : (0);
    for (int j = 0; j < ny_; ++j) {
      int jn = (j > 0) ? (j - 1) : (ny_ - 1);
      int jp = (j < ny_ - 1) ? (j + 1) : (0);

      f[i][j][0] = fcol[i][j][0];
      f[ip][j][1] = fcol[i][j][1];
      f[i][jp][2] = fcol[i][j][2];
      f[in][j][3] = fcol[i][j][3];
      f[i][jn][4] = fcol[i][j][4];
      f[ip][jp][5] = fcol[i][j][5];
      f[in][jp][6] = fcol[i][j][6];
      f[in][jn][7] = fcol[i][j][7];
      f[ip][jn][8] = fcol[i][j][8];
    }
  }
}

lbgrid::~lbgrid() {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) delete[] f[i][j];
    delete[] f[i];
  }
  delete[] f;
}
