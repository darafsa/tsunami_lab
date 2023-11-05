/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Unit tests for the one-dimensional wave propagation patch.
 **/
#include <catch2/catch.hpp>
#include "WavePropagation1d.h"

TEST_CASE( "Test the 1d wave propagation solver.", "[WaveProp1d]" ) {
  /*
   * Test case:
   *
   *   Single dam break problem between cell 49 and 50.
   *     left | right
   *       10 | 8
   *        0 | 0
   *
   *   Elsewhere steady state.
   *
   * The net-updates at the respective edge are given as
   * (see derivation in Roe solver):
   *    left          | right
   *     9.394671362  | -9.394671362
   *    -88.25985     | -88.25985
   */

  // construct solver and setup a dambreak problem
  tsunami_lab::patches::WavePropagation1d m_waveProp( 100 );

  for( std::size_t l_ce = 0; l_ce < 50; l_ce++ ) {
    m_waveProp.setHeight( l_ce,
                          0,
                          10 );
    m_waveProp.setMomentumX( l_ce,
                             0,
                             0 );
  }
  for( std::size_t l_ce = 50; l_ce < 100; l_ce++ ) {
    m_waveProp.setHeight( l_ce,
                          0,
                          8 );
    m_waveProp.setMomentumX( l_ce,
                             0,
                             0 );
  }

  // set outflow boundary condition
  m_waveProp.setGhostOutflow();

  // perform a time step
  m_waveProp.timeStep( 0.1, tsunami_lab::patches::WavePropagation::FWave );

  // steady state
  for( std::size_t l_ce = 0; l_ce < 49; l_ce++ ) {
    REQUIRE( m_waveProp.getHeight()[l_ce]   == Approx(10) );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }

  // dam-break
  REQUIRE( m_waveProp.getHeight()[49]   == Approx(10 - 0.1 * 9.394671362) );
  REQUIRE( m_waveProp.getMomentumX()[49] == Approx( 0 + 0.1 * 88.25985) );

  REQUIRE( m_waveProp.getHeight()[50]   == Approx(8 + 0.1 * 9.394671362) );
  REQUIRE( m_waveProp.getMomentumX()[50] == Approx(0 + 0.1 * 88.25985) );

  // steady state
  for( std::size_t l_ce = 51; l_ce < 100; l_ce++ ) {
    REQUIRE( m_waveProp.getHeight()[l_ce]   == Approx(8) );
    REQUIRE( m_waveProp.getMomentumX()[l_ce] == Approx(0) );
  }
}

TEST_CASE("Test the 1d wave propagation FWave solver (Shock-Shock Problem).", "[WaveProp1dFWaveShockShock]")
{
  /*
   * @brief test state from middle_states.csv (Shock-Shock Problem)
   *
   * h_l = 3042.136044684769
   * h_r = 3042.136044684769
   * hu_l = -27.52440428024561
   * hu_r = 27.52440428024561
   * h* = 3041.976718035753
   */

  // construct solver and setup a Shock-Shock problem
  tsunami_lab::patches::WavePropagation1d m_waveProp(100);

  for (std::size_t l_ce = 0; l_ce < 50; l_ce++)
  {
    m_waveProp.setHeight( l_ce,
                          0,
                          3042.136044684769);
    m_waveProp.setMomentumX( l_ce,
                             0,
                             -27.52440428024561);
  }
  for (std::size_t l_ce = 50; l_ce < 100; l_ce++)
  {
    m_waveProp.setHeight( l_ce,
                          0,
                          3042.136044684769);
    m_waveProp.setMomentumX( l_ce,
                             0,
                             27.52440428024561);
  }

  // set outflow boundary condition
  m_waveProp.setGhostOutflow();

  // perform a time step
  m_waveProp.timeStep( 0.1, tsunami_lab::patches::WavePropagation::FWave );

  // test for h*
  REQUIRE( m_waveProp.getHeight()[50]   == Approx(3041.976718035753) );
}

TEST_CASE("Test the 1d wave propagation FWave solver (Rare-Rare Problem", "[WaveProp1dFWaveRareRare]")
{
  /**
   * @brief test state from middle_states.csv (Rare-Rare Problem)
   *
   * h_l = 7589.71304876485
   * h_r = 7589.71304876485
   * hu_l = -138.9853242339589
   * hu_r = 138.9853242339589
   * h* = 7589.203700916305
   */

  // construct solver and setup a Rare-Rare problem
  tsunami_lab::patches::WavePropagation1d m_waveProp(100);

  for (std::size_t l_ce = 0; l_ce < 50; l_ce++)
  {
    m_waveProp.setHeight(l_ce,
                         0,
                         7589.71304876485);
    m_waveProp.setMomentumX(l_ce,
                            0,
                            -138.9853242339589);
  }
  for (std::size_t l_ce = 50; l_ce < 100; l_ce++)
  {
    m_waveProp.setHeight(l_ce,
                         0,
                         7589.71304876485);
    m_waveProp.setMomentumX(l_ce,
                            0,
                            138.9853242339589);
  }

  // set outflow boundary condition
  m_waveProp.setGhostOutflow();

  // perform a time step
  m_waveProp.timeStep( 0.1, tsunami_lab::patches::WavePropagation::FWave );

  // test for h*
  REQUIRE(m_waveProp.getHeight()[50] == Approx(7589.203700916305));
}