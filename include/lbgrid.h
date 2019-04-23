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
  void initialize_density(double rho);

  //! Returns the density function at a node in a certain direction
  //! \param[in] i number of node in the x direction
  //! \param[in] j number of node in the y direction
  //! \param[in] k direction (0 to 8)
  double density_function(int i, int j, int k) const { return f[i][j][k]; };

  //! Returns the total mass
  double sum_density();

  //! Destructor
  ~lbgrid();

 private:
  int nx_;
  int ny_;
  const int Q = 9;
  double*** f = new double**[nx_];
};

#endif
