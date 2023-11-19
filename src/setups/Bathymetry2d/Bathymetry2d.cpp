/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * One-dimensional dam break problem.
 **/
#include "Bathymetry2d.h"

tsunami_lab::setups::Bathymetry2d::Bathymetry2d( real in_heightInner, real in_heightOuter, real in_radiusDam, real in_xMax, real in_yMax, real in_scaling ) {
	heightInner = in_heightInner;
	heightOuter = in_heightOuter;
	radiusDam = in_radiusDam * in_scaling;
	scaling = in_scaling;

	centerDam[0] = (in_xMax/2) * scaling;
	centerDam[1] = (in_yMax/2) * scaling;
}

tsunami_lab::t_real tsunami_lab::setups::Bathymetry2d::getHeight( t_real in_x, t_real in_y ) const {
	real distanceFromCenter = sqrt((centerDam[0]-in_x)*(centerDam[0]-in_x) + (centerDam[1]-in_y)*(centerDam[1]-in_y));
	// t_real distanceFromCenter = sqrt(in_x*in_x + in_y*in_y);
	if( distanceFromCenter < radiusDam ) {
		return heightInner;
	}
	else {
		return heightOuter;
	}
}

tsunami_lab::t_real tsunami_lab::setups::Bathymetry2d::getMomentumX( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::Bathymetry2d::getMomentumY( t_real, t_real ) const {
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::Bathymetry2d::getBathymetry( t_real in_x, t_real in_y ) const {
	if ( in_x > 7.5 && in_x < 8 && in_y > 7.5 && in_y < 8 ) {
		return -1;
	}
	return -2;
}