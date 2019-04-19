#ifndef LBGRID_H
#define LBGRID_H
#include <iostream>

class lbgrid {
 public:
  lbgrid(int nx, int ny);
  unsigned nx() const { return nx_; };
  unsigned ny() const { return ny_; };
  void initialize_density(double rho);
  double density_function(int i, int j, int k) const { return f[i][j][k]; };
  double sum_density();
  ~lbgrid();

 private:
  int nx_;
  int ny_;
  const int Q = 9;
  double*** f = new double**[nx_];
 };

#endif
