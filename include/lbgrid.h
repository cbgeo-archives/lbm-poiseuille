#ifndef LBGRID_H
#define LBGRID_H
#include <cmath>
#include <iostream>

//! \brief A class that stores the information about the grid
//! \details This class stores number of nodes in each direction and the density
//! functions at each node.
class lbgrid {
 public:
  //! \param[in] nx number of nodes in x direction
  //! \param[in] ny number of nodes in y direction
  lbgrid(int nx, int ny);

  //! Returns nx of the grid
  unsigned nx() const { return nx_; };

  //! Returns ny of the grid
  unsigned ny() const { return ny_; };

  //! Initializes the density functions
  //! \param[in] rho density
  void initialize_density(double rho_0);

  //! Returns the density function at a node in a certain direction
  //! \param[in] i node number in the x direction
  //! \param[in] j node number in the y direction
  //! \param[in] k lattice velocity direction (0 to 8)
  double f_(int i, int j, int k) const { return f[i][j][k]; };

  //! Returns the total mass
  double sum_density();

  //! Computes the macroscopic variable including velocity and density
  //! \param[in] Fx body force in the x direction
  void compute_macro_var(double Fx, double Fy);

  //! Return the macroscopic velocity in the x direction at a node
  //! \param[in] i node number in the x direction
  //! \param[in] j node number in the y direction
  double ux_(int i, int j) { return ux[i][j]; };

  //! Returns the macroscopic velocity in the y direction at a node
  //! \param[in] i node number in the x direction
  //! \param[in] j node number in the y direction
  double uy_(int i, int j) { return uy[i][j]; };

  //! Return the macroscopic velocities in the x direction at a cross section
  //! \param[in] i node number in the x direction
  double* uxsection(int i) { return ux[i]; };

  //! Computes the equilibrium density function
  void equilibrium_density();

  //! Returns the equilibrium density function at a node in certain direction
  //! \param[in] i node number in the x direction
  //! \param[in] j node number in the y direction
  //! \param[in] k lattice velocity direction (0 to 8)
  double feq_(int i, int j, int k) { return feq[i][j][k]; };

  //! Computes the density function after collision
  //! \param[in] tau normalized relaxation time
  //! \param[in] Fx body force in the x direction
  //! \param[in] Fy body force in the y direction
  void collision(double tau, double Fx, double Fy);

  //! Return the density function after collision at a node in a certain
  //! direction \param[in] i node number in the x direction \param[in] j node
  //! number in the y direction \param[in] k lattice velocity direction (0 to 8)
  double fcol_(int i, int j, int k) { return fcol[i][j][k]; };

  //! Streams the density functions
  void streaming();

  //! Calculates the analytical solution
  //! \param[in] tau normalized relaxation time
  //! \param[in] nt number of time step
  void analytical_solution(double tau, int nt, double Fx);

  //! Returns the analytical velocity in the x direction at a cross section
  //! \param[in] i node number in the x direction
  double* u_an_section(int i) { return u_an[i]; };

  //! Writes VTK
  //! \param[in] nt number of time step
  void write_vtk(int nt);

  //! Destructor
  ~lbgrid();

 private:
  int nx_;
  int ny_;
  const int Q = 9;
  const int ex[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
  const int ey[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
  const double w[9] = {4. / 9.,  1. / 9.,  1. / 9.,  1. / 9., 1. / 9.,
                       1. / 36., 1. / 36., 1. / 36., 1. / 36.};
  const double w0 = 4. / 9.;
  const double w1_4 = 1. / 9.;
  const double w5_8 = 1. / 36.;
  double rho;
  double*** f = new double**[nx_];
  double** ux = new double*[nx_];
  double** uy = new double*[nx_];
  double*** feq = new double**[nx_];
  double*** fcol = new double**[nx_];
  double** u_an = new double*[nx_];
};
#endif
