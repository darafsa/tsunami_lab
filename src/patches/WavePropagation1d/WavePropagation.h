/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Base class of the wave propagation patches.
 **/
#ifndef TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION
#define TSUNAMI_LAB_PATCHES_WAVE_PROPAGATION

#include "../../constants.h"

namespace tsunami_lab {
  namespace patches {
    class WavePropagation;
  }
}

class tsunami_lab::patches::WavePropagation {
  public:
    /**
     * @brief Virtual destructor for base class.
     **/
    virtual ~WavePropagation(){};

    /**
     * @brief Performs a time step.
     *
     * @param in_scaling scaling of the time step.
     **/
    virtual void timeStep( real in_scaling, Solver in_solver ) = 0; 

	 /**
	  * @brief Sets the values of the ghost cells according to outflow boundary conditions.
	  * 
	  * @param in_boundary boundary type to use (outflow/reflective); 0: boundary left side, 1: boundary right side.
	  */
    virtual void setGhostOutflow(Boundary in_boundary[2]) = 0;

    /**
     * @brief Gets the stride in y-direction. x-direction is stride-1.
     *
     * @return stride in y-direction.
     **/
    virtual idx getStride() = 0;

    /**
     * @brief Gets cells' water heights.
     *
     * @return water heights.
     */
    virtual real const * getHeight() = 0;

    /**
     * @brief Gets the cells' momenta in x-direction.
     *
     * @return momenta in x-direction.
     **/
    virtual real const * getMomentumX() = 0;

    /**
     * @brief Gets the cells' momenta in y-direction.
     *
     * @return momenta in y-direction.
     **/
    virtual real const * getMomentumY() = 0;

	 /**
     * @brief Gets the cells' bathymetry.
     *
     * @return bathymetry.
     **/
    virtual real const * getBathymetry() = 0;

    /**
     * @brief Sets the height of the cell to the given value.
     *
     * @param in_x id of the cell in x-direction.
     * @param in_y id of the cell in y-direction.
     * @param in_height water height.
     **/
    virtual void setHeight( idx  in_x,
                            idx  in_y,
                            real in_height ) = 0;

    /**
     * @brief Sets the momentum in x-direction to the given value.
     *
     * @param in_x id of the cell in x-direction.
     * @param in_y id of the cell in y-direction.
     * @param in_momentum momentum in x-direction.
     **/
    virtual void setMomentumX( idx  in_x,
                               idx  in_y,
                               real in_momentum ) = 0;

    /**
     * @brief Sets the momentum in y-direction to the given value.
     *
     * @param in_x id of the cell in x-direction.
     * @param in_y id of the cell in y-direction.
     * @param in_momentumVertical momentum in y-direction.
     **/
    virtual void setMomentumY( idx  in_x,
                               idx  in_y,
                               real in_momentumVertical ) = 0;

	 /**
	  * @brief Sets the bathymetry to the given value.
	  * 
	  * @param in_x id of the cell in x-direction.
     * @param in_y id of the cell in y-direction.
     * @param in_bathymetry bathymetry.
	  **/
	 virtual void setBathymetry( idx in_x,
	 									  idx in_y,
										  real in_bathymetry ) = 0;
};

#endif