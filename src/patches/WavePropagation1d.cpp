/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional wave propagation patch.
 **/
#include "WavePropagation1d.h"
#include "../solvers/FWave.h"
#include "../solvers/Roe.h"

using namespace tsunami_lab::patches;

WavePropagation1d::WavePropagation1d( idx in_cellCount ) {
  cellCount = in_cellCount;

  // allocate memory including a single ghost cell on each side
  for( unsigned short step = 0; step < 2; step++ ) {
    height[step] = new real[ cellCount + 2 ];
    momentum[step] = new real[ cellCount + 2 ];
  }
  bathymetry = new real[ cellCount + 2 ];

  // init to zero
  for( unsigned short step = 0; step < 2; step++ ) {
    for( idx cell = 0; cell < cellCount; cell++ ) {
      height[step][cell] = 0;
      momentum[step][cell] = 0;
		if(step == 0) {
			bathymetry[cell] = 0;
		}
    }
  }
}

WavePropagation1d::~WavePropagation1d() {
  for( unsigned short i = 0; i < 2; i++ ) {
    delete[] height[i];
    delete[] momentum[i];
  }
  delete[] bathymetry;
}

void WavePropagation1d::timeStep( real in_scaling, Solver in_solver ) {
  // pointers to old and new data
  real * heightOld = height[step];
  real * momentumOld = momentum[step];

  step = (step+1) % 2;
  real * heightNew =  height[step];
  real * momentumNew = momentum[step];

  // init new cell quantities
  for( idx cell = 1; cell < cellCount+1; cell++ ) {
    heightNew[cell] = heightOld[cell];
    momentumNew[cell] = momentumOld[cell];
  }

  // iterate over edges and update with Riemann solutions
  for( idx edge = 0; edge < cellCount+1; edge++ ) {
    // determine left and right cell-id
    idx cellLeft = edge;
    idx cellRight = edge+1;

    // compute net-updates
    real netUpdates[2][2];
	 
	 real stateLeft[3] = { heightOld[cellLeft], momentumOld[cellLeft], bathymetry[cellLeft] };
	 real stateRight[3] = { heightOld[cellRight], momentumOld[cellRight], bathymetry[cellRight] };

	 if ( in_solver == FWAVE ) {
		solvers::FWave::netUpdates( stateLeft, 
	 										 stateRight, 
                              	 netUpdates[0],
                              	 netUpdates[1] );
	 } else {
		solvers::Roe::netUpdates( stateLeft[0], 
	 									  stateRight[0], 
	 									  stateLeft[1], 
	 									  stateRight[1], 
                                netUpdates[0],
                                netUpdates[1] );
	 }

    // update the cells' quantities
    heightNew[cellLeft]  -= in_scaling * netUpdates[0][0];
    momentumNew[cellLeft] -= in_scaling * netUpdates[0][1];

    heightNew[cellRight]  -= in_scaling * netUpdates[1][0];
    momentumNew[cellRight] -= in_scaling * netUpdates[1][1];
  }
}

void WavePropagation1d::setGhostOutflow( Boundary boundary[2] ) {
  real * heightLocal = height[step];
  real * momentumLocal = momentum[step];
  real * bathymetryLocal = bathymetry;

  // set left boundary
  if(boundary[0] == OUTFLOW) {
	 heightLocal[0] = heightLocal[1];
	 momentumLocal[0] = momentumLocal[1];
	 bathymetryLocal[0] = bathymetryLocal[1];
  } else if (boundary[1] == REFLECTING) {
	 heightLocal[0] = 0;
	 momentumLocal[0] = 0;
	 bathymetryLocal[0] = heightLocal[1]+1;
  }

  // set right boundary
  if(boundary[1] == OUTFLOW) {
	 heightLocal[cellCount+1] = heightLocal[cellCount];
	 momentumLocal[cellCount+1] = momentumLocal[cellCount];
	 bathymetryLocal[cellCount+1] = bathymetryLocal[cellCount];
  } else if(boundary[1] == REFLECTING) {
	 heightLocal[cellCount+1] = 0;
	 momentumLocal[cellCount+1] = 0;
	 bathymetryLocal[cellCount+1] = bathymetryLocal[cellCount]+1;
  }
}