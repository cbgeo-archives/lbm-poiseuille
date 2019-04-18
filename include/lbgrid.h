#ifndef LBGRID_H
#define LBGRID_H
#include <iostream>

class lbgrid {
 private:
  int nx_;
  int ny_;
  const int Q = 9;
  double*** f = new double**[nx_];
  double total_mass = 0.0;

 public:
  lbgrid(int nx, int ny);
  void initialize_density(double rho);
  double density_function(int i, int j, int k);
  double sum_density();
  ~lbgrid();
};

#endif
