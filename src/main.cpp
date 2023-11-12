/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Entry-point for simulations.
 **/
#include "io/Csv.h"
#include "patches/WavePropagation1d.h"
#include "setups/DamBreak1d/DamBreak1d.h"
#include "setups/RareRare1d/RareRare1d.h"
#include "setups/ShockShock1d/ShockShock1d.h"
#include "setups/Bathymetry1d/Bathymetry1d.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

int main(int in_argc, char *in_argv[]) {
  // number of cells in x- and y-direction
  tsunami_lab::idx xCount = 0;
  tsunami_lab::idx yCount = 1;

  // set cell size
  tsunami_lab::real cellSize = 1;

  // solver type
  tsunami_lab::patches::WavePropagation1d::Solver solverType;

  std::cout << "####################################" << std::endl;
  std::cout << "### Tsunami Lab                  ###" << std::endl;
  std::cout << "###                              ###" << std::endl;
  std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
  std::cout << "####################################" << std::endl;

  if (in_argc < 4) {
    std::cerr << "invalid number of arguments, usage:" << std::endl;
    std::cerr << "  ./build/tsunami_lab CELLS SOLVER SETUP height velocity" << std::endl;
    std::cerr << "where CELLS is the number of cells in x-direction, "
                 "SOLVER the solver type [FWAVE, ROE] and SETUP the setup to "
                 "use [DAMBREAK, RARE, SHOCK]."
              << std::endl;
    return EXIT_FAILURE;
  } else {
    xCount = atoi(in_argv[1]);
    if (xCount < 1) {
      std::cerr << "invalid number of cells" << std::endl;
      return EXIT_FAILURE;
    }
    cellSize = 10.0 / xCount;

    std::string solverArg = in_argv[2];
    if (solverArg == "FWAVE") {
      solverType = tsunami_lab::patches::WavePropagation::FWave;
    } else if (solverArg == "ROE") {
      solverType = tsunami_lab::patches::WavePropagation::Roe;
    } else {
      std::cerr << "invalid solver type. Please use either ROE or FWAVE"
                << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::cout << "runtime configuration" << std::endl;
  std::cout << "  number of cells in x-direction: " << xCount << std::endl;
  std::cout << "  number of cells in y-direction: " << yCount << std::endl;
  std::cout << "  cell size:                      " << cellSize << std::endl;

  // construct setup
  std::string setupArg = in_argv[3];
  tsunami_lab::setups::Setup *setup;
  tsunami_lab::real height = 10;
  tsunami_lab::real momentum = 50;
  if (setupArg == "DAMBREAK") {
    setup = new tsunami_lab::setups::DamBreak1d(10, 5, 5);
  } else if (setupArg == "RARE") {
	 if (in_argc > 4) {
		height = std::stof(in_argv[4]);
		momentum = std::stof(in_argv[5]) * height;
	 }
    setup = new tsunami_lab::setups::RareRare1d(height, momentum, 5);
  } else if (setupArg == "SHOCK") {
	 if (in_argc > 4) {
		height = std::stof(in_argv[4]);
		momentum = std::stof(in_argv[5]) * height;
	 }
    setup = new tsunami_lab::setups::ShockShock1d(height, momentum, 5);
  } else if(setupArg == "BATHYMETRY") {
	 setup = new tsunami_lab::setups::Bathymetry1d(10, 5, 5);
  } else {
    std::cerr << "invalid setup type. Please use either DAMBREAK, RARE or SHOCK"
              << std::endl;
    return EXIT_FAILURE;
  }

  // construct solver
  tsunami_lab::patches::WavePropagation *waveProp;
  waveProp = new tsunami_lab::patches::WavePropagation1d(xCount);

  // maximum observed height in the setup
  tsunami_lab::real heightMax =
      std::numeric_limits<tsunami_lab::real>::lowest();

  // set up solver
  for (tsunami_lab::idx cellY = 0; cellY < yCount; cellY++) {
    tsunami_lab::real y = cellY * cellSize;

    for (tsunami_lab::idx cellX = 0; cellX < xCount; cellX++) {
      tsunami_lab::real x = cellX * cellSize;

      // get initial values of the setup
      tsunami_lab::real height = setup->getHeight(x, y);
      heightMax = std::max(height, heightMax);

      tsunami_lab::real momentumX = setup->getMomentumX(x, y);
      tsunami_lab::real momentumY = setup->getMomentumY(x, y);
      tsunami_lab::real bathymetry = setup->getBathymetry(x, y);

      // set initial values in wave propagation solver
      waveProp->setHeight(cellX, cellY, height);

      waveProp->setMomentumX(cellX, cellY, momentumX);

      waveProp->setMomentumY(cellX, cellY, momentumY);

      waveProp->setBathymetry(cellX, cellY, bathymetry);
    }
  }

  // derive maximum wave speed in setup; the momentum is ignored
  tsunami_lab::real speedMax = std::sqrt(9.81 * heightMax);

  // derive constant time step; changes at simulation time are ignored
  tsunami_lab::real dt = 0.5 * cellSize / speedMax;

  // derive scaling for a time step
  tsunami_lab::real scaling = dt / cellSize;

  // set up time and print control
  tsunami_lab::idx timeStep = 0;
  tsunami_lab::idx nOut = 0;
  tsunami_lab::real endTime = 3;
  tsunami_lab::real simTime = 0;

  std::cout << "entering time loop" << std::endl;

  // iterate over time
  while (simTime < endTime) {
    if (timeStep % 25 == 0) {
      std::cout << "  simulation time / #time steps: " << simTime << " / "
                << timeStep << std::endl;

      std::string path = "solution_" + std::to_string(nOut) + ".csv";
      std::cout << "  writing wave field to " << path << std::endl;

      std::ofstream file;
      file.open(path);

      tsunami_lab::io::Csv::write(cellSize, 
											 xCount, 
											 1, 1, 
											 waveProp->getHeight(),
                                  waveProp->getBathymetry(), 
											 waveProp->getMomentumX(), 
											 nullptr, 
											 file);
      file.close();
      nOut++;
    }

	 tsunami_lab::patches::WavePropagation::Boundary boundary[2] = {
		tsunami_lab::patches::WavePropagation::Outflow, 
		tsunami_lab::patches::WavePropagation::Outflow };
    waveProp->setGhostOutflow(boundary);
    waveProp->timeStep(scaling, solverType);

    timeStep++;
    simTime += dt;
  }

  std::cout << "finished time loop" << std::endl;

  // free memory
  std::cout << "freeing memory" << std::endl;
  delete setup;
  delete waveProp;

  std::cout << "finished, exiting" << std::endl;
  return EXIT_SUCCESS;
}
