/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 *
 * 1d shock shock problem setup.
 *
 **/

#include "ShockShock1d.h"

tsunami_lab::setups::ShockShock1d::ShockShock1d(t_real i_height,
                                                t_real i_momentum,
                                                t_real i_midPos)
{
    m_height = i_height;
    m_momentum = i_momentum;
    m_midPos = i_midPos;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getHeight(t_real,
                                                                 t_real) const
{
  return m_height;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getMomentumX(t_real i_x,
                                                                    t_real) const
{
  if (i_x < m_midPos)
    {
        return m_momentum;
    }
  else
    {
        return -m_momentum;
    }
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getMomentumY(t_real,
                                                                    t_real) const
{
  return 0;
}

tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getBathymetry(t_real,
                                                                  	t_real) const
{
  return 0;
}