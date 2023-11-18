/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional subcritical case.
 **/
#ifndef TSUNAMI_LAB_SETUPS_SUBCRITICAL_1D_H
#define TSUNAMI_LAB_SETUPS_SUBCRITICAL_1D_H

#include "../Setup.h"

namespace tsunami_lab {
  namespace setups {
    class Subcritical1d;
  }
}

/**
 * @brief 1d subcritical case setup.
 **/
class tsunami_lab::setups::Subcritical1d: public Setup {
  

  public:
    /**
     * @brief Gets the water height at a given point.
     *
     * @param in_x x-coordinate of the queried point.
     * @return height at the given point.
     **/
    t_real getHeight( t_real in_x,
                      t_real      ) const;

    /**
     * @brief Gets the momentum in x-direction.
     *
     * @param in_x x-coordinate of the queried point.
     * @return momentum in x-direction.
     **/
    t_real getMomentumX( t_real in_x,
                         t_real     ) const;

    /**
     * @brief Gets the momentum in y-direction.
     *
     * @return momentum in y-direction.
     **/
    t_real getMomentumY( t_real,
                         t_real ) const;

	 /**
     * @brief Gets the bathymetry at a given point.
     *
     * @param in_x x-coordinate of the queried point.
     * @return bathymetry.
     **/
    t_real getBathymetry( t_real in_x,
                          t_real     ) const;

};

#endif