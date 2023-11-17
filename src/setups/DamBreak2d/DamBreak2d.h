/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Two-dimensional dam break problem.
 **/
#ifndef TSUNAMI_LAB_SETUPS_DAM_BREAK_2D_H
#define TSUNAMI_LAB_SETUPS_DAM_BREAK_2D_H

#include <cmath>
#include "../Setup.h"

namespace tsunami_lab {
	namespace setups {
		class DamBreak2d;
	}
}

/**
 * @brief 2d dam break setup.
 **/
class tsunami_lab::setups::DamBreak2d: public Setup {
	private:
		//! height on the inner side
		real heightInner = 0;

		//! height on the outer side
		real heightOuter = 0;

		//! radius of the dam
		real radiusDam = 0;

		//! center point of dam; 0: x, 1: y
		idx centerDam[2];

	public:
		/**
		 * @brief Constructor.
		 *
		 * @param in_heightInner water height on the inner side of the dam.
		 * @param in_heightOuter water height on the outer side of the dam.
		 * @param in_radiusDam radius of the dam.
		 * @param in_cellCountX cell count in x-direction.
		 * @param in_cellCountY cell count in y-direction.
		**/
		DamBreak2d( real in_heightInner, real in_heightOuter, real in_radiusDam, idx in_cellCountX, idx in_cellCountY  );

		/**
		 * @brief Gets the water height at a given point.
		 *
		 * @param in_x x-coordinate of the queried point.
		 * @return height at the given point.
		**/
		real getHeight( t_real in_x, t_real y ) const;

		/**
		 * @brief Gets the momentum in x-direction.
		 *
		 * @return momentum in x-direction.
		**/
		real getMomentumX( t_real, t_real ) const;

		/**
		 * @brief Gets the momentum in y-direction.
		 *
		 * @return momentum in y-direction.
		**/
		real getMomentumY( t_real, t_real ) const;

		/**
		 * @brief Gets the bathymetry at a given point.
		 *
		 * @return bathymetry.
		**/
   	real getBathymetry( t_real, t_real ) const;

};

#endif