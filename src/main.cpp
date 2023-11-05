/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Entry-point for simulations.
 **/
#include "io/Csv.h"
#include "patches/WavePropagation1d.h"
#include "setups/DamBreak1d.h"
#include "setups/RareRare1d.h"
#include "setups/ShockShock1d.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

int main(int i_argc, char *i_argv[]) {
  // number of cells in x- and y-direction
  tsunami_lab::t_idx l_nx = 0;
  tsunami_lab::t_idx l_ny = 1;

  // set cell size
  tsunami_lab::t_real l_dxy = 1;

  // solver type
  tsunami_lab::patches::WavePropagation1d::Solver l_solverType;

  std::cout << "####################################" << std::endl;
  std::cout << "### Tsunami Lab                  ###" << std::endl;
  std::cout << "###                              ###" << std::endl;
  std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
  std::cout << "####################################" << std::endl;

  if (i_argc < 4) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    std::cerr << "  ./build/tsunami_lab CELLS SOLVER SETUP height velocity" << std::endl;
    std::cerr << "where CELLS is the number of cells in x-direction, "
                 "SOLVER the solver type [FWAVE, ROE] and SETUP the setup to "
                 "use [DAMBREAK, RARE, SHOCK]."
              << std::endl;
    return EXIT_FAILURE;
  } else {
    l_nx = atoi(i_argv[1]);
    if (l_nx < 1) {
      std::cerr << "invalid number of cells" << std::endl;
      return EXIT_FAILURE;
    }
    l_dxy = 10.0 / l_nx;

    std::string l_solver = i_argv[2];
    if (l_solver == "FWAVE") {
      l_solverType = tsunami_lab::patches::WavePropagation::FWave;
    } else if (l_solver == "ROE") {
      l_solverType = tsunami_lab::patches::WavePropagation::Roe;
    } else {
      std::cerr << "invalid solver type. Please use either ROE or FWAVE"
                << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << l_nx << std::endl;
  std::cout << "  number of cells in y-direction: " << l_ny << std::endl;
  std::cout << "  cell size:                      " << l_dxy << std::endl;

  // construct setup
  std::string l_setupArg = i_argv[3];
  tsunami_lab::setups::Setup *l_setup;
  tsunami_lab::t_real l_height = 10;
  tsunami_lab::t_real l_momentum = 50;
  if (l_setupArg == "DAMBREAK") {
    l_setup = new tsunami_lab::setups::DamBreak1d(10, 5, 5);
  } else if (l_setupArg == "RARE") {
	 if (i_argc > 4) {
		l_height = std::stof(i_argv[4]);
		l_momentum = std::stof(i_argv[5]) * l_height;
	 }
    l_setup = new tsunami_lab::setups::RareRare1d(l_height, l_momentum, 5);
  } else if (l_setupArg == "SHOCK") {
	 if (i_argc > 4) {
		l_height = std::stof(i_argv[4]);
		l_momentum = std::stof(i_argv[5]) * l_height;
	 }
    l_setup = new tsunami_lab::setups::ShockShock1d(l_height, l_momentum, 5);
  } else {
    std::cerr << "invalid setup type. Please use either DAMBREAK, RARE or SHOCK"
              << std::endl;
    return EXIT_FAILURE;
  }

  // construct solver
  tsunami_lab::patches::WavePropagation *l_waveProp;
  l_waveProp = new tsunami_lab::patches::WavePropagation1d(l_nx);

  // maximum observed height in the setup
  tsunami_lab::t_real l_hMax =
      std::numeric_limits<tsunami_lab::t_real>::lowest();

  // set up solver
  for (tsunami_lab::t_idx l_cy = 0; l_cy < l_ny; l_cy++) {
    tsunami_lab::t_real l_y = l_cy * l_dxy;

    for (tsunami_lab::t_idx l_cx = 0; l_cx < l_nx; l_cx++) {
      tsunami_lab::t_real l_x = l_cx * l_dxy;

      // get initial values of the setup
      tsunami_lab::t_real l_h = l_setup->getHeight(l_x, l_y);
      l_hMax = std::max(l_h, l_hMax);

      tsunami_lab::t_real l_hu = l_setup->getMomentumX(l_x, l_y);
      tsunami_lab::t_real l_hv = l_setup->getMomentumY(l_x, l_y);

      // set initial values in wave propagation solver
      l_waveProp->setHeight(l_cx, l_cy, l_h);

      l_waveProp->setMomentumX(l_cx, l_cy, l_hu);

      l_waveProp->setMomentumY(l_cx, l_cy, l_hv);
    }
  }

  // derive maximum wave speed in setup; the momentum is ignored
  tsunami_lab::t_real l_speedMax = std::sqrt(9.81 * l_hMax);

  // derive constant time step; changes at simulation time are ignored
  tsunami_lab::t_real l_dt = 0.5 * l_dxy / l_speedMax;

  // derive scaling for a time step
  tsunami_lab::t_real l_scaling = l_dt / l_dxy;

  // set up time and print control
  tsunami_lab::t_idx l_timeStep = 0;
  tsunami_lab::t_idx l_nOut = 0;
  tsunami_lab::t_real l_endTime = 1.25;
  tsunami_lab::t_real l_simTime = 0;

  std::cout << "entering time loop" << std::endl;

  // iterate over time
  while (l_simTime < l_endTime) {
    if (l_timeStep % 25 == 0) {
      std::cout << "  simulation time / #time steps: " << l_simTime << " / "
                << l_timeStep << std::endl;

      std::string l_path = "solution_" + std::to_string(l_nOut) + ".csv";
      std::cout << "  writing wave field to " << l_path << std::endl;

      std::ofstream l_file;
      l_file.open(l_path);

      tsunami_lab::io::Csv::write(l_dxy, l_nx, 1, 1, l_waveProp->getHeight(),
                                  l_waveProp->getMomentumX(), nullptr, l_file);
      l_file.close();
      l_nOut++;
    }

    l_waveProp->setGhostOutflow();
    l_waveProp->timeStep(l_scaling, l_solverType);

    l_timeStep++;
    l_simTime += l_dt;
  }

  std::cout << "finished time loop" << std::endl;

  // free memory
  std::cout << "freeing memory" << std::endl;
  delete l_setup;
  delete l_waveProp;

  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
