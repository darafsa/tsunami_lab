/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional shock shock problem setup.
 **/
#ifndef TSUNAMI_LAB_SETUPS_SHOCK_SHOCK_REFLECTIVE_1D_H
#define TSUNAMI_LAB_SETUPS_SHOCK_SHOCK_REFLECTIVE_1D_H

#include "../Setup.h"

namespace tsunami_lab {
  namespace setups {
    class ShockShockReflective1d;
  }
}

/**
 * 1d dam break setup.
 **/
class tsunami_lab::setups::ShockShockReflective1d: public Setup {
  private:
    //! height of the water 
    t_real m_height = 0;
    
    //! momentum of the water
    t_real m_momentum = 0;

    //! location of the middle point
    t_real m_midPos = 0;

	 Boundary m_boundaryLeft = OUTFLOW;
	 Boundary m_boundaryRight = OUTFLOW;

  public:
    /**
     * Constructor.
     *
     * @param i_height water height wich is the same on both sides.
     * @param i_momentum water momentum which is the same on both sides.
     * @param i_midPos the middle position where both waves collide/separate.
     **/
    ShockShockReflective1d( t_real i_height,
               				 t_real i_momentum,
               				 t_real i_midPos );

    /**
     * Gets the water height at a given point.
     *
     * @param i_x x-coordinate of the queried point.
     * @return height at the given point.
     **/
    t_real getHeight( t_real i_x,
                      t_real      ) const;

    /**
     * Gets the momentum in x-direction.
     *
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real,
                         t_real ) const;

    /**
     * Gets the momentum in y-direction.
     *
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real,
                         t_real ) const;
							
	 /**
     * @brief Gets the bathymetry at a given point.
     *
     * @param i_x x-coordinate of the queried point.
     * @return bathymetry.
     **/
    t_real getBathymetry( t_real i_x,
                          t_real ) const;
};

#endif