/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 *
 * 1d supercritical case setup.
 *
 **/

#include "Supercritical1d.h"
#include <cmath>

tsunami_lab::t_real tsunami_lab::setups::Supercritical1d::getHeight(t_real in_x, t_real) const {
	if (in_x >= 0 && in_x <= 25) {
		return -getBathymetry(in_x, 0);
	}
	else {
		return 0;
	}
}

tsunami_lab::t_real tsunami_lab::setups::Supercritical1d::getMomentumX(t_real in_x, t_real) const {
	if (in_x >= 0 && in_x <= 25) {
		return 0.18;
	}
	else {
		return 0;
	}
}

tsunami_lab::t_real tsunami_lab::setups::Supercritical1d::getMomentumY(t_real, t_real) const {
	return 0;
}

tsunami_lab::t_real tsunami_lab::setups::Supercritical1d::getBathymetry(t_real in_x, t_real) const {
	if (in_x > 8 && in_x < 12) {
		return (-0.13 - 0.05 * (in_x - 10) * (in_x - 10));
	}
	else { 
		return -0.33;
	}
}