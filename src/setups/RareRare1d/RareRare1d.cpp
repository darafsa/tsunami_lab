/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 *
 * One-dimensional rare rare problem setup.
 *
 **/

#include "RareRare1d.h"

tsunami_lab::setups::RareRare1d::RareRare1d(t_real i_height,
                                            t_real i_momentum,
                                            t_real i_midPos)
{
    m_height = i_height;
    m_momentum = i_momentum;
    m_midPos = i_midPos;
}

tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getHeight(t_real,
                                                               t_real) const
{
  return m_height;
}

tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getMomentumX(t_real i_x,
                                                                  t_real) const
{
  if (i_x < m_midPos)
    {
        return -m_momentum;
    }
  else
    {
        return m_momentum;
    }
}

tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getMomentumY(t_real,
                                                                  t_real) const
{
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getBathymetry(t_real, 
																						 t_real) const {
	return 0;
}