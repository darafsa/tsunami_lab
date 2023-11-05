/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional rare rare problem setup.
 **/
#ifndef TSUNAMI_LAB_SETUPS_RARE_RARE_1D_H
#define TSUNAMI_LAB_SETUPS_RARE_RARE_1D_H

#include "Setup.h"

namespace tsunami_lab {
  namespace setups {
    class RareRare1d;
  }
}

/**
 * 1d dam break setup.
 **/
class tsunami_lab::setups::RareRare1d: public Setup {
  private:
    //! height of the water 
    t_real m_height = 0;
    
    //! momentum of the water
    t_real m_momentum = 0;

    //! location of the middle point
    t_real m_midPos = 0;

  public:
    /**
     * Constructor.
     *
     * @param i_heightLeft water height on the left side of the dam.
     * @param i_heightRight water height on the right side of the dam.
     * @param i_locationDam location (x-coordinate) of the dam.
     **/
    RareRare1d( t_real i_heightLeft,
                t_real i_heightRight,
                t_real i_locationDam );

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

};

#endif