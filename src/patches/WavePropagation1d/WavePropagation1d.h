/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional wave propagation patch.
 **/
#ifndef TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_1D
#define TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION_1D

#include "../WavePropagation.h"

namespace tsunami_lab {
  namespace patches {
    class WavePropagation1d;
  }
}

class tsunami_lab::patches::WavePropagation1d: public WavePropagation {
  private:
    //! current step which indicates the active values in the arrays below
    unsigned short step = 0;

    //! number of cells discretizing the computational domain
    idx cellCount = 0;

    //! water heights for the current and next time step for all cells
    real * height[2] = { nullptr, nullptr };

    //! momenta for the current and next time step for all cells
    real * momentum[2] = { nullptr, nullptr };

    //! bathymetry for all cells
    real * bathymetry = nullptr;

	 //! minmal bathymetry depth
	 real dy = -20;

  public:
    /**
     * @brief Constructs the 1d wave propagation solver.
     *
     * @param in_cellCount number of cells.
     **/
    WavePropagation1d( idx in_cellCount );

    /**
     * @brief Destructor which frees all allocated memory.
     **/
    ~WavePropagation1d();

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
	  * @param in_boundary boundary type to use (outflow/reflective); 0: boundary left side, 1: boundary right side.
	  */
    void setGhostOutflow(Boundary in_boundary[2]);

    /**
     * @brief Gets the stride in y-direction. x-direction is stride-1.
     *
     * @return stride in y-direction.
     **/
    idx getStride(){
      return cellCount+2;
    }

    /**
     * @brief Gets cells' water heights.
     *
     * @return water heights.
     */
    real const * getHeight(){
      return height[step]+1;
    }

    /**
     * @brief Gets the cells' momenta in x-direction.
     *
     * @return momenta in x-direction.
     **/
    real const * getMomentumX(){
      return momentum[step]+1;
    }

    /**
     * @brief Dummy function which returns a nullptr.
     **/
    real const * getMomentumY(){
      return nullptr;
    }

	 real const * getBathymetry(){
		return bathymetry+1;
	 }

    /**
     * @brief Sets the height of the cell to the given value.
     *
     * @param in_x id of the cell in x-direction.
     * @param in_height water height.
     **/
    void setHeight( idx  in_x,
	 					  idx,
                    real in_height ) {
      height[step][in_x+1] = in_height;
    }

    /**
     * @brief Sets the momentum in x-direction to the given value.
     *
     * @param in_x id of the cell in x-direction.
     * @param in_momentumHorizontal momentum in x-direction.
     **/
    void setMomentumX( idx  in_x,
	 						  idx,
                       real in_momentumHorizontal ) {
      momentum[step][in_x+1] = in_momentumHorizontal;
    }

    /**
     * @brief Dummy function since there is no y-momentum in the 1d solver.
     **/
    void setMomentumY( idx,
                       idx,
                       real ) {};
	
	 /**
	  * @brief Sets the bathymetry to the given value.
	  * 
	  * @param in_x id of the cell in x-direction.
     * @param in_bathymetry bathymetry.
	  **/
	 void setBathymetry( idx in_x,
	 							idx,
	 							real in_bathymetry ) {
		// if(in_bathymetry > dy && in_bathymetry < 0) {
		// 	in_bathymetry = dy;
		// }
		bathymetry[in_x + 1] = in_bathymetry;
	};
};

#endif