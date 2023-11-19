/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 *
 * 1d subcritical case setup.
 *
 **/

#include "Subcritical1d.h"

tsunami_lab::t_real tsunami_lab::setups::Subcritical1d::getHeight(t_real in_x, t_real) const {
	if (in_x >= 0 && in_x <= 25) {
		return -getBathymetry(in_x, 0);
	}
	else {
		return 0;
	}
}

tsunami_lab::t_real tsunami_lab::setups::Subcritical1d::getMomentumX(t_real in_x, t_real) const {
	if (in_x >= 0 && in_x <= 25) {
		return 4.42;
	}
	else {
		return 0;
	}
}

tsunami_lab::t_real tsunami_lab::setups::Subcritical1d::getMomentumY(t_real, t_real) const {
	return 0;
}

tsunami_lab::t_real tsunami_lab::setups::Subcritical1d::getBathymetry(t_real in_x, t_real) const {
	if (in_x > 8 && in_x < 12) {
		return (-1.8 - 0.05 * (in_x - 10) * (in_x - 10));
	}
	else {
		return -2;
	}
}