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
      feq[i][j] = new double[Q]();
      fcol[i][j] = new double[Q]();
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

void lbgrid::compute_macro_var(double dt, double Fx, double Fy) {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      rho = f[i][j][0] + f[i][j][1] + f[i][j][2] + f[i][j][3] + f[i][j][4] +
            f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];
      std::cout<<"rho is = "<<rho<<std::endl;
      ux[i][j] = (f[i][j][1] + f[i][j][5] + f[i][j][8] -
                  (f[i][j][3] + f[i][j][6] + f[i][j][7])) /
                 rho;
      uy[i][j] = (f[i][j][2] + f[i][j][5] + f[i][j][6] -
                  (f[i][j][4] + f[i][j][7] + f[i][j][8])) /
                 rho;
    }
  }
}

void lbgrid::equilibrium_density(double dt) {

  double a1 = 3.;
  double a2 = 9. / 2.;
  double a3 = 3. / 2.;

  for (int i = 0; i < nx_; ++i) {
    for (int j = 1; j < ny_-1; ++j) {

      const double uxsq = ux[i][j] * ux[i][j];
      const double uysq = uy[i][j] * uy[i][j];
      const double usq = uxsq + uysq;

      rho = f[i][j][0] + f[i][j][1] + f[i][j][2] + f[i][j][3] + f[i][j][4] +
            f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];

      feq[i][j][0] = w0 * rho * (1. - a3 * usq);
      feq[i][j][1] = w1_4 * rho * (1. + a1 * ux[i][j] + a2 * uxsq - a3 * usq);
      feq[i][j][2] = w1_4 * rho * (1. + a1 * uy[i][j] + a2 * uysq - a3 * usq);
      feq[i][j][3] = w1_4 * rho * (1. - a1 * ux[i][j] + a2 * uxsq - a3 * usq);
      feq[i][j][4] = w1_4 * rho * (1. - a1 * uy[i][j] + a2 * uysq - a3 * usq);
      feq[i][j][5] =
          w5_8 * rho *
          (1. + a1 * (ux[i][j] + uy[i][j]) + a2 * (uxsq + uysq) - a3 * usq);
      feq[i][j][6] =
          w5_8 * rho *
          (1. + a1 * (-ux[i][j] + uy[i][j]) + a2 * (-uxsq + uysq) - a3 * usq);
      feq[i][j][7] =
          w5_8 * rho *
          (1. + a1 * (-ux[i][j] - uy[i][j]) + a2 * (-uxsq - uysq) - a3 * usq);
      feq[i][j][8] =
          w5_8 * rho *
          (1. + a1 * (ux[i][j] - uy[i][j]) + a2 * (uxsq - uysq) - a3 * usq);
    }
  }
}

void lbgrid::collision(double tau, double dt, double Fx, double Fy) {
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {

      double S[8];
      S[0] = 0.;
      S[1] = (1 - 0.5 * dt / tau) * w1_4 *(3 * Fx);
      S[2]=(1-0.5*dt/tau)*w1_4*(3*Fy);
      S[3]=(1-0.5*dt/tau)*w1_4*(-3*Fx);
      S[4]=(1-0.5*dt/tau)*w1_4*(-3*Fy);
      S[5] =(1-0.5*dt/tau)*w5_8*(3*Fx+3*Fy); 
      S[6] =(1-0.5*dt/tau)*w5_8*(-3*Fx+3*Fy); 
      S[7] =(1-0.5*dt/tau)*w5_8*(-3*Fx-3*Fy); 
      S[8] =(1-0.5*dt/tau)*w5_8*(3*Fx-3*Fy); 

      for (int q = 0; q < Q; ++q) {
        fcol[i][j][q] =
            f[i][j][q] - (f[i][j][q] - feq[i][j][q]) / tau * dt + S[q] * dt;
      }
    }
  }
}

void lbgrid::streaming() {
  for (int i = 0; i < nx_; ++i) {
    f[i][0][2]=fcol[i][0][4];
    f[i][0][5]=fcol[i][0][7];
    f[i][0][6]=fcol[i][0][8];
    f[i][ny_-1][4]=fcol[i][ny_-1][2];
    f[i][ny_-1][7]=fcol[i][ny_-1][5];
    f[i][ny_-1][8]=fcol[i][ny_-1][6];

    int in = (i > 0) ? (i - 1) : (nx_ - 1);
    int ip = (i < nx_ - 1) ? (i + 1) : (0);
    for (int j = 0; j < ny_; ++j) {
      int jn = j-1;
      int jp = j + 1;

      f[i][j][0] = fcol[i][j][0];
      f[ip][j][1] = fcol[i][j][1];
      f[in][j][3] = fcol[i][j][3];

      if (j!=ny_-1) {
      f[i][jp][2] = fcol[i][j][2];
      f[ip][jp][5] = fcol[i][j][5];
      f[in][jp][6] = fcol[i][j][6];
      }

      if (j!=0) {
      f[i][jn][4] = fcol[i][j][4];
      f[in][jn][7] = fcol[i][j][7];
      f[ip][jn][8] = fcol[i][j][8];
      }
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
