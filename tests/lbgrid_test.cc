#include "catch.hpp"

#include "lbgrid.h"

TEST_CASE("Functionality of the lbgrid class is checked", "[lbgrid]") {

  const int nx = 2, ny = 3;
  const double rho = 5.;
  const double Tolerance = 1.E-7;

  std::unique_ptr<lbgrid> grid = std::make_unique<lbgrid>(nx, ny);
  //  REQUIRE(grid->nx==2);
  //  REQUIRE(grid->ny==3);
  REQUIRE(grid->density_function(0, 0, 0) == Approx(0.).epsilon(Tolerance));

  grid->initialize_density(rho);
  REQUIRE(grid->density_function(0, 0, 0) ==
          Approx(4. * rho / 9.).epsilon(Tolerance));

  REQUIRE(grid->sum_density() == Approx(nx * ny * rho).epsilon(Tolerance));
}
