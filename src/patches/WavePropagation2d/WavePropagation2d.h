/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Two-dimensional wave propagation patch.
 **/
#ifndef TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D
#define TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_2D

#include "../WavePropagation.h"

namespace tsunami_lab {
	namespace patches {
		class WavePropagation2d;
	}
}

class tsunami_lab::patches::WavePropagation2d: public WavePropagation {
	private:
		//! current step which indicates the active values in the arrays below
		unsigned short step = 0;

		//! number of cells discretizing the computational domain
		idx cellCountX = 0;
		idx cellCountY = 0;

		//! water heights for the current and next time step for all cells
		real ** height[2] = { nullptr, nullptr };

		//! momenta for the current and next time step for all cells in x-direction
		real ** momentumX[2] = { nullptr, nullptr };

		//! momenta for the current and next time step for all cells in y-direction
		real ** momentumY[2] = { nullptr, nullptr };

		//! array used to convert 2d array into 1d
		real * array1d;

		//! bathymetry for all cells
		real ** bathymetry;

		real const * linearizeArray(real ** array2d) {
			for(idx y = 0; y < cellCountY; y++) {
				for(idx x = 0; x < cellCountX; x++) {
					array1d[x + y*cellCountX] = array2d[x+1][y+1];
				}
			}
			return array1d;
		};

		void copyGhostCellsOutflow( real ** out_grid );
		void copyGhostCellsReflecting( real ** out_grid, real in_value );

	public:
		/**
		 * @brief Constructs the 2d wave propagation solver.
		 *
		 * @param in_cellCountX number of cells in x-direction.
		 * @param in_cellCountY number of cells in y-direction.
		 **/
		WavePropagation2d( idx in_cellCountX, idx in_cellCountY );

		/**
		 * @brief Destructor which frees all allocated memory.
		 **/
		~WavePropagation2d();

		/**
		 * @brief Performs a time step.
		 *
		 * @param in_scaling scaling of the time step (dt / dx).
		 * @param in_solver solver type to use (Roe / FWave)
		 **/
		void timeStep( real in_scaling, Solver in_solver );

		/**
		* @brief Sets the values of the ghost cells according to outflow boundary conditions.
		* 
		* @param in_boundary boundary type to use (outflow/reflective); 0: boundary -x, 1: boundary x, 2: boundary -y, 3: boundary y.
		*/
		void setGhostOutflow(Boundary in_boundary[4]);

		/**
		 * @brief Gets the stride in y-direction. x-direction is stride-1.
		 *
		 * @return stride in y-direction.
		 **/
		idx getStride(){
			return 0;
		}

		/**
		 * @brief Gets cells' water heights.
		 *
		 * @return water heights.
		 */
		real const * getHeight(){
			return linearizeArray(height[step]);
		}

		/**
		 * @brief Gets the cells' momenta in x-direction.
		 *
		 * @return momenta in x-direction.
		 **/
		real const * getMomentumX(){
			return linearizeArray(momentumX[step]);
		}

		/**
		 * @brief Gets the cells' momenta in y-direction.
		 *
		 * @return momenta in y-direction.
		 **/
		real const * getMomentumY(){
			return linearizeArray(momentumY[step]);
		}

		/**
		 * @brief Get the bathymetry.
		 * 
		 * @return bathymetry.
		 */
		real const * getBathymetry(){
			return linearizeArray(bathymetry);
		}

		/**
		 * @brief Sets the height of the cell to the given value.
		 *
		 * @param in_x id of the cell in x-direction.
		 * @param in_x id of the cell in y-direction.
		 * @param in_height water height.
		 **/
		void setHeight( idx in_x, idx in_y, real in_height ) {
			height[step][in_x+1][in_y+1] = in_height;
		}

		/**
		 * @brief Sets the momentum in x-direction to the given value.
		 *
		 * @param in_x id of the cell in x-direction.
		 * @param in_y id of the cell in y-direction.
		 * @param in_momentumX momentum in x-direction.
		 **/
		void setMomentumX( idx in_x, idx in_y, real in_momentumX ) {
			momentumX[step][in_x+1][in_y+1] = in_momentumX;
		}

		/**
		 * @brief Sets the momentum in y-direction to the given value.
		 *
		 * @param in_x id of the cell in x-direction.
		 * @param in_y id of the cell in y-direction.
		 * @param in_momentumY momentum in y-direction.
		 **/
		void setMomentumY( idx in_x, idx in_y, real in_momentumY ) {
			momentumY[step][in_x+1][in_y+1] = in_momentumY;
		}
	
		/**
		 * @brief Sets the bathymetry to the given value.
		 * 
		 * @param in_x id of the cell in x-direction.
		 * @param in_y id of the cell in y-direction.
		 * @param in_bathymetry bathymetry.
		**/
		void setBathymetry( idx in_x, idx in_y, real in_bathymetry ) {
			bathymetry[in_x+1][in_y+1] = in_bathymetry;
		};
};

#endif