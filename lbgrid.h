#ifndef LBGRID_H
#define LBGRID_H
#include <iostream>

class lbgrid {
private:
  int nx_;
  int ny_;
  const int Q = 9;
  double ***f = new double **[nx_];

public:
  lbgrid(int nx, int ny);
  void initialize_density(double rho);
  void print_f();
};

#endif
