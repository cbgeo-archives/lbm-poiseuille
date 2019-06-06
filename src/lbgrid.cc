#include "lbgrid.h"
#include <iomanip>
#include <omp.h>

lbgrid::lbgrid(int nx, int ny) : nx_{nx}, ny_{ny} {
  for (unsigned i = 0; i < nx_; ++i) {
    f[i] = new double*[ny_]();
    ux[i] = new double[ny_]();
    uy[i] = new double[ny_]();
    feq[i] = new double*[ny_]();
    fcol[i] = new double*[ny_]();
    u_an[i] = new double[ny_]();
    for (unsigned j = 0; j < ny_; ++j) {
      f[i][j] = new double[Q]();
      feq[i][j] = new double[Q]();
      fcol[i][j] = new double[Q]();
    }
  }
}

void lbgrid::initialize_density(double rho_0) {
  int i, j, k;
#pragma omp parallel for default(shared) private(i, j, k) schedule(dynamic)
  for (i = 0; i < nx_; ++i) {
    for (j = 0; j < ny_; ++j) {
      for (k = 0; k < Q; ++k) {
        f[i][j][k] = w[k] * rho_0;
      }
    }
  }
}

double lbgrid::sum_density() {
  double total_mass = 0.0;
  int i, j, k;
#pragma omp parallel for default(shared) private(i,j,k) schedule(dynamic) reduction(+:total_mass)
  for (i = 0; i < nx_; ++i) {
    for (j = 0; j < ny_; ++j) {
      for (k = 0; k < Q; ++k) {
        total_mass += f[i][j][k];
      }
    }
  }
  return total_mass;
}

void lbgrid::compute_macro_var(double Fx, double Fy) {
  double fy;
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      fy = 0.;
      if (j == ny_ / 2) fy = Fy;
      rho = ux[i][j] = uy[i][j] = 0.;
      for (int k = 0; k < Q; ++k) {
        rho += f[i][j][k];
        ux[i][j] += f[i][j][k] * ex[k];
        uy[i][j] += f[i][j][k] * ey[k];
      }
      ux[i][j] = ux[i][j] / rho + Fx / 2;
      uy[i][j] = uy[i][j] / rho + fy / 2;
    }
  }
}

void lbgrid::equilibrium_density() {

  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      rho = 0.;
      for (int k = 0; k < Q; ++k) rho += f[i][j][k];
      for (int k = 0; k < Q; ++k) {
        feq[i][j][k] = w[k] * rho *
                       (1. + 3. * (ux[i][j] * ex[k] + uy[i][j] * ey[k]) +
                        9. / 2. * (ux[i][j] * ex[k] + uy[i][j] * ey[k]) *
                            (ux[i][j] * ex[k] + uy[i][j] * ey[k]) -
                        3. / 2. * (ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j]));
      }
    }
  }
}

void lbgrid::collision(double tau, double Fx, double Fy) {
  double fy;
  double S[9];
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      fy = 0.;
      if (j == ny_ / 2) fy = Fy;
      for (int k = 0; k < Q; ++k) {
        S[k] = (1. - 0.5 / tau) * w[k] *
               (3 * (Fx * ex[k] + fy * ey[k]) +
                9 * (Fx * ex[k] + fy * ey[k]) *
                    (ex[k] * ux[i][j] + ey[k] * uy[i][j]) -
                3 * (ux[i][j] * Fx + uy[i][j] * fy));

        fcol[i][j][k] = f[i][j][k] - (f[i][j][k] - feq[i][j][k]) / tau + S[k];
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
      f[in][j][3] = fcol[i][j][3];
      if (j != ny_ - 1) {
        f[i][jp][2] = fcol[i][j][2];
        f[ip][jp][5] = fcol[i][j][5];
        f[in][jp][6] = fcol[i][j][6];
      }

      if (j != 0) {
        f[i][jn][4] = fcol[i][j][4];
        f[in][jn][7] = fcol[i][j][7];
        f[ip][jn][8] = fcol[i][j][8];
      }
    }
  }

  for (int i = 0; i < nx_; ++i) {
    f[i][0][2] = fcol[i][0][4];
    f[i][0][5] = fcol[i][0][7];
    f[i][0][6] = fcol[i][0][8];
    f[i][ny_ - 1][4] = fcol[i][ny_ - 1][2];
    f[i][ny_ - 1][7] = fcol[i][ny_ - 1][5];
    f[i][ny_ - 1][8] = fcol[i][ny_ - 1][6];
  }
}

/*void lbgrid::analytical_solution(double tau, int nt, double Fx) {
  double y, t;
  double a;
  double nu = (2 * tau - 1) / 6;
  double umax = Fx * nx_ * nx_ / 8 / nu;
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      y = -1. + (2. * j + 1.) / ny_;
      t = (nt * (tau - 0.5)) / (3. * (0.5 * ny_) * (0.5 * ny_));
      a = 0.;
      for (int n = 0; n < 50; ++n) {
        a += pow(-1, n) * cos(M_PI * y * (n + 0.5)) / pow(M_PI, 3) /
             pow((n + 0.5), 3) * exp(-pow(M_PI, 2) * t * pow((n + 0.5), 2));
      }
      u_an[i][j] = umax * ((1. - y * y) - 4. * a);
    }
  }
}*/

void lbgrid::write_vtk(int nt) {
  int i, j;
  double pasxyz;
  double P, u_x, u_y, uan;
  char filename[25];
  FILE* sortie;

  sprintf(filename, "lbm%.6i.vtk", nt);
  pasxyz = 1.;
  sortie = fopen(filename, "w");
  fprintf(sortie, "# vtk DataFile Version 2.0\n");
  fprintf(sortie, "Sortie domaine LB+LINK t %e\n", nt);
  fprintf(sortie, "ASCII\n");
  fprintf(sortie, "DATASET RECTILINEAR_GRID\n");
  fprintf(sortie, "DIMENSIONS %d %d 1\n", nx_, ny_);
  fprintf(sortie, "X_COORDINATES %d float\n", nx_);
  for (i = 0; i < nx_; ++i) fprintf(sortie, "%e ", (float)i * pasxyz);
  fprintf(sortie, "\n");
  fprintf(sortie, "Y_COORDINATES %d float\n", ny_);
  for (i = 0; i < ny_; ++i) fprintf(sortie, "%e ", (float)i * pasxyz);
  fprintf(sortie, "\n");
  fprintf(sortie, "Z_COORDINATES 1 float\n");
  fprintf(sortie, "0\n");

  // For LB
  fprintf(sortie, "POINT_DATA %d\n", nx_ * ny_);
  fprintf(sortie, "VECTORS VecVelocity float\n");

  for (j = 0; j < ny_; ++j) {
    for (i = 0; i < nx_; ++i) {
      u_x = ux[i][j];
      u_y = uy[i][j];
      fprintf(sortie, "%.8lf %.8lf 0.\n", u_x, u_y);
    }
  }

  /* fprintf(sortie, "VECTORS AnalVelocity float\n");

   for (j = 0; j < ny_; ++j) {
     for (i = 0; i < nx_; ++i) {
       uan = u_an[i][j];
       fprintf(sortie, "%.4lf 0. 0.\n", uan);
     }
   }*/
  fclose(sortie);
}

lbgrid::~lbgrid() {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      delete[] f[i][j];
      delete[] fcol[i][j];
      delete[] feq[i][j];
    }
    delete[] f[i];
    delete[] fcol[i];
    delete[] feq[i];
    delete[] ux[i];
    delete[] uy[i];
    delete[] u_an[i];
  }
  delete[] f;
  delete[] fcol;
  delete[] feq;
  delete[] ux;
  delete[] uy;
  delete[] u_an;
}
