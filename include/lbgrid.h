#ifndef LBGRID_H
#define LBGRID_H
#include <iostream>

//! \brief A class that stores the information about the grid
//! \details This class stores number of nodes in each direction and the density
//! functions at each node.
class lbgrid {
 public:
  //! Constructs a grid with assigned density functions at each node
  //! \param[in] nx number of nodes in x direction
  //! \param[in] ny number of nodes in y direction
  lbgrid(int nx, int ny);

  //! Returns nx of the grid
  unsigned nx() const { return nx_; };

  //! Returns ny of the grid
  unsigned ny() const { return ny_; };

  //! Calculates the initial density functions
  //! \param[in] rho density
  void initialize_density(double rho_0);

  //! Returns the density function at a node in a certain direction
  //! \param[in] i number of node in the x direction
  //! \param[in] j number of node in the y direction
  //! \param[in] k direction (0 to 8)
  double f_(int i, int j, int k) const { return f[i][j][k]; };

  //! Returns the total mass
  double sum_density();

  void compute_macro_var (double dt, double Fx, double Fy);

  double ux_(int i, int j) {return ux[i][j];};
  double uy_(int i, int j) {return uy[i][j];}; 
  double * uxsection(int i) {return ux[i];};

  void equilibrium_density(double dt);
  double feq_(int i, int j , int k) {return feq[i][j][k];};

  void collision(double tau, double dt, double Fx, double Fy);
  double fcol_(int i, int j, int k) {return fcol[i][j][k];};
  void streaming();

  //! Destructor
  ~lbgrid();

 private:
  int nx_;
  int ny_;
  const int Q = 9;
  const double w0 = 4./9.;
  const double w1_4 =1./9.;
  const double w5_8 = 1./36.;
  double rho;
  double*** f = new double**[nx_];
  double** ux=new double*[nx_];
  double** uy=new double*[nx_];
  double*** feq= new double**[nx_];
  double*** fcol= new double**[nx_];
};
#endif
