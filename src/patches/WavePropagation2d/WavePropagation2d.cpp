/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional wave propagation patch.
 **/
#include "WavePropagation2d.h"
#include "../../solvers/FWave.h"
#include "../../solvers/Roe.h"

using namespace tsunami_lab::patches;

WavePropagation2d::WavePropagation2d( idx in_cellCountX, idx in_cellCountY ) {
	cellCountX = in_cellCountX;
	cellCountY = in_cellCountY;

	// allocate memory including a single ghost cell on each side
	for( unsigned short step = 0; step < 2; step++ ) {
		height[step] = new real*[ cellCountX + 2 ];
		momentumX[step] = new real*[ cellCountX + 2 ];
		momentumY[step] = new real*[ cellCountX + 2 ];
		for ( idx x = 0; x < cellCountX + 2; x++) {
			height[step][x] = new real[ cellCountY + 2];
			momentumX[step][x] = new real[ cellCountY + 2];
			momentumY[step][x] = new real[ cellCountY + 2];
		}
	}
	bathymetry = new real*[ cellCountX + 2 ];
	for ( idx x = 0; x < cellCountX + 2; x++) {
		bathymetry[x] = new real[ cellCountY + 2];
	}
	array1d = new real[ cellCountX * cellCountY ];

	// init to zero
	for( unsigned short step = 0; step < 2; step++ ) {
		for( idx y = 0; y < cellCountY + 2; y++ ) {
			for( idx x = 0; x < cellCountX + 2; x++ ) {
				height[step][x][y] = 0;
				momentumX[step][x][y] = 0;
				momentumY[step][x][y] = 0;
				if( step == 0 ) {
					bathymetry[x][y] = 0;
				}
			}
		}
	}
	for ( idx cell = 0; cell < cellCountX * cellCountY; cell++ ) {
		array1d[cell] = 0;
	}
}

WavePropagation2d::~WavePropagation2d() {
	for( unsigned short step = 0; step < 2; step++ ) {
		for ( idx x = 0; x < cellCountX + 2; x++) {
			delete height[step][x];
			delete momentumX[step][x];
			delete momentumY[step][x];
		}
		delete[] height[step];
		delete[] momentumX[step];
		delete[] momentumY[step];
	}
	for ( idx x = 0; x < cellCountX + 2; x++) {
		delete bathymetry[x];
	}
	delete[] bathymetry;
}

void WavePropagation2d::timeStep( real in_scaling, Solver in_solver ) {
	// pointers to old and new data
	real ** heightOld = height[step];
	real ** momentumXOld = momentumX[step];
	real ** momentumYOld = momentumY[step];

	step = (step+1) % 2;
	real ** heightNew =	height[step];
	real ** momentumXNew = momentumX[step];
	real ** momentumYNew = momentumY[step];

	// init new cell quantities
	for( idx y = 1; y < cellCountY + 1; y++ ) {
		for( idx x = 1; x < cellCountX + 1; x++) {
			heightNew[x][y] = heightOld[x][y];
			momentumXNew[x][y] = momentumXOld[x][y];
			momentumYNew[x][y] = momentumYOld[x][y];
		}
	}

	// iterate over edges and update with Riemann solutions in x-direction
	for( idx y = 0; y < cellCountY + 2; y++ ) {
		for( idx edgeX = 0; edgeX < cellCountX + 1; edgeX++ ) {
			// determine cell-id
			idx cellLeft = edgeX;
			idx cellRight = edgeX+1;

			// compute net-updates
			real netUpdates[2][2];
		
			real stateLeft[3] = { heightOld[cellLeft][y], momentumXOld[cellLeft][y], bathymetry[cellLeft][y] };
			real stateRight[3] = { heightOld[cellRight][y], momentumXOld[cellRight][y], bathymetry[cellRight][y] };

			if(bathymetry[cellLeft][y] > 0) {
				stateLeft[0] = stateRight[0];
				stateLeft[1] = -stateRight[1];
				stateLeft[2] = stateRight[2];
			}

			if(bathymetry[cellRight][y] > 0) {
				stateRight[0] = stateLeft[0];
				stateRight[1] = -stateLeft[1];
				stateRight[2] = stateLeft[2];
			}

			if ( in_solver == FWAVE ) {
				solvers::FWave::netUpdates( stateLeft, stateRight, netUpdates[0], netUpdates[1] );
			} else {
				solvers::Roe::netUpdates( stateLeft[0], stateRight[0], stateLeft[1], stateRight[1], netUpdates[0], netUpdates[1] );
			}

			// update the cells' quantities
			heightNew[cellLeft][y] -= in_scaling * netUpdates[0][0];
			momentumXNew[cellLeft][y] -= in_scaling * netUpdates[0][1];

			heightNew[cellRight][y]	-= in_scaling * netUpdates[1][0];
			momentumXNew[cellRight][y] -= in_scaling * netUpdates[1][1];
		}
	}

	// iterate over edges and update with Riemann solutions in y-direction
	for( idx x = 0; x < cellCountX + 2; x++ ) {
		for( idx edgeY = 0; edgeY < cellCountX + 1; edgeY++ ) {
			// determine cell-id
			idx cellTop = edgeY+1;
			idx cellBottom = edgeY;

			// compute net-updates
			real netUpdates[2][2];
		
			real stateLeft[3] = { heightOld[x][cellBottom], momentumYOld[x][cellBottom], bathymetry[x][cellBottom] };
			real stateRight[3] = { heightOld[x][cellTop], momentumYOld[x][cellTop], bathymetry[x][cellTop] };

			if(bathymetry[x][cellBottom] > 0) {
				stateLeft[0] = stateRight[0];
				stateLeft[1] = -stateRight[1];
				stateLeft[2] = stateRight[2];
			}

			if(bathymetry[x][cellTop] > 0) {
				stateRight[0] = stateLeft[0];
				stateRight[1] = -stateLeft[1];
				stateRight[2] = stateLeft[2];
			}

			if ( in_solver == FWAVE ) {
				solvers::FWave::netUpdates( stateLeft, stateRight, netUpdates[0], netUpdates[1] );
			} else {
				solvers::Roe::netUpdates( stateLeft[0], stateRight[0], stateLeft[1], stateRight[1], netUpdates[0], netUpdates[1] );
			}

			// update the cells' quantities
			heightNew[x][cellBottom] -= in_scaling * netUpdates[0][0];
			momentumYNew[x][cellBottom] -= in_scaling * netUpdates[0][1];

			heightNew[x][cellTop] -= in_scaling * netUpdates[1][0];
			momentumYNew[x][cellTop] -= in_scaling * netUpdates[1][1];
		}
	}
}

void WavePropagation2d::copyGhostCellsOutflow( real ** out_grid ) {
	idx xMax = cellCountX+1;
	idx yMax = cellCountY+1;

	for( idx x = 1; x < xMax; x++ ) {
		out_grid[x][0] = out_grid[x][1];
		out_grid[x][yMax] = out_grid[x][yMax-1];
	}

	for( idx y = 1; y < yMax; y++ ) {
		out_grid[0][y] = out_grid[1][y];
		out_grid[xMax][y] = out_grid[xMax-1][y];
	}

	out_grid[0][0] = out_grid[1][1];
	out_grid[xMax][0] = out_grid[xMax-1][1];
	out_grid[0][yMax] = out_grid[1][yMax-1];
	out_grid[xMax][yMax] = out_grid[xMax-1][yMax-1];
}

void WavePropagation2d::copyGhostCellsReflecting( real ** out_grid, real in_value ) {
	idx xMax = cellCountX+1;
	idx yMax = cellCountY+1;

	for( idx x = 1; x < xMax; x++ ) {
		out_grid[x][0] = in_value;
		out_grid[x][yMax] = in_value;
	}

	for( idx y = 1; y < yMax; y++ ) {
		out_grid[0][y] = in_value;
		out_grid[xMax][y] = in_value;
	}

	out_grid[0][0] = in_value;
	out_grid[xMax][0] = in_value;
	out_grid[0][yMax] = in_value;
	out_grid[xMax][yMax] = in_value;
}

void WavePropagation2d::setGhostOutflow( Boundary in_boundary[2] ) {
	// set left boundary
	if(in_boundary[0] == OUTFLOW) {
		copyGhostCellsOutflow( height[step] );
		copyGhostCellsOutflow( momentumX[step] );
		copyGhostCellsOutflow( momentumY[step] );
		copyGhostCellsOutflow( bathymetry );
	} else if (in_boundary[0] == REFLECTING) {
		copyGhostCellsReflecting( height[step], 0 );
		copyGhostCellsReflecting( momentumX[step], 0 );
		copyGhostCellsReflecting( momentumY[step], 0 );
		copyGhostCellsReflecting( bathymetry, 20 );
	}
}