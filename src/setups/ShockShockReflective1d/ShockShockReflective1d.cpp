/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 *
 * 1d shock shock problem setup.
 *
 **/

#include "ShockShockReflective1d.h"

tsunami_lab::setups::ShockShockReflective1d::ShockShockReflective1d( t_real i_height,
                                                							t_real i_momentum,
                                                							t_real i_midPos )
{
    m_height = i_height;
    m_momentum = i_momentum;
    m_midPos = i_midPos;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShockReflective1d::getHeight( t_real,
                                                                		 		 t_real ) const
{
  return m_height;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShockReflective1d::getMomentumX( t_real,
                                                                    				 t_real ) const
{
	return m_momentum;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShockReflective1d::getMomentumY( t_real,
                                                                     			 t_real ) const
{
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShockReflective1d::getBathymetry( t_real,
                                                                  				  t_real) const
{
  return 0;
}