/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional dam break problem.
 **/
#include "Bathymetry1d.h"
#include <math.h>

tsunami_lab::setups::Bathymetry1d::Bathymetry1d( t_real i_heightLeft,
                                             	 t_real i_heightRight,
                                             	 t_real i_locationDam ) {
  m_heightLeft = i_heightLeft;
  m_heightRight = i_heightRight;
  m_locationDam = i_locationDam;
}

tsunami_lab::t_real tsunami_lab::setups::Bathymetry1d::getHeight( t_real i_x,
                                                               	t_real ) const {
  if( i_x < m_locationDam ) {
    return m_heightLeft - getBathymetry(i_x, 0);
  }
  else {
    return m_heightRight - getBathymetry(i_x, 0);
  }
}

tsunami_lab::t_real tsunami_lab::setups::Bathymetry1d::getMomentumX( t_real,
                                                                   	t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::Bathymetry1d::getMomentumY( t_real,
                                                                   	t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::Bathymetry1d::getBathymetry( t_real i_x, 
																						  	 t_real ) const {
	// return -(sin(i_x) + t_real(1))-1;

	if (i_x > 7.5 && i_x < 7.6) {
		return -1;
	}
	return -2;
}