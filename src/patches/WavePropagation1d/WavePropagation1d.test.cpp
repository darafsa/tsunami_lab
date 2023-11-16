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
  tsunami_lab::patches::WavePropagation1d waveProp( 100 );

  for( std::size_t cell = 0; cell < 50; cell++ ) {
    waveProp.setHeight( cell,
                        0,
                        10 );
    waveProp.setMomentumX( cell,
                           0,
                           0 );
	 waveProp.setBathymetry( cell,
	 								 0,
									 0 );
  }
  for( std::size_t cell = 50; cell < 100; cell++ ) {
    waveProp.setHeight( cell,
                        0,
                        8 );
    waveProp.setMomentumX( cell,
                           0,
                           0 );
	 waveProp.setBathymetry( cell,
	 								 0,
									 0 );
  }

  // set outflow boundary condition
  tsunami_lab::Boundary boundary[2] = {
		tsunami_lab::OUTFLOW, 
		tsunami_lab::OUTFLOW };
  waveProp.setGhostOutflow( boundary );

  // perform a time step
  waveProp.timeStep( 0.1, tsunami_lab::FWAVE );

  // steady state
  for( std::size_t cell = 0; cell < 49; cell++ ) {
    REQUIRE( waveProp.getHeight()[cell]   == Approx(10) );
    REQUIRE( waveProp.getMomentumX()[cell] == Approx(0) );
  }

  // dam-break
  REQUIRE( waveProp.getHeight()[49]   == Approx(10 - 0.1 * 9.394671362) );
  REQUIRE( waveProp.getMomentumX()[49] == Approx( 0 + 0.1 * 88.25985) );

  REQUIRE( waveProp.getHeight()[50]   == Approx(8 + 0.1 * 9.394671362) );
  REQUIRE( waveProp.getMomentumX()[50] == Approx(0 + 0.1 * 88.25985) );

  // steady state
  for( std::size_t cell = 51; cell < 100; cell++ ) {
    REQUIRE( waveProp.getHeight()[cell]   == Approx(8) );
    REQUIRE( waveProp.getMomentumX()[cell] == Approx(0) );
  }
}

TEST_CASE("Test the 1d wave propagation FWave solver (Shock-Shock Problem).", "[WaveProp1dFWaveShockShock]")
{
  /*
   * @brief test state from middle_states.csv (Shock-Shock Problem)
   *
   * h_l = 8899.326826472694
   * h_r = 8899.326826472694
   * hu_l = 122.0337839252433
   * hu_r = -122.0337839252433
   * h* = 8899.739847378269
   */

  // construct solver and setup a Shock-Shock problem
  tsunami_lab::patches::WavePropagation1d waveProp( 100 );

  for (std::size_t cell = 0; cell < 50; cell++)
  {
    waveProp.setHeight( cell,
                        0,
                        8899.326826472694 );
    waveProp.setMomentumX( cell,
                           0,
                           122.0337839252433 );
	 waveProp.setBathymetry( cell,
	 								 0,
									 0 );
  }
  for (std::size_t cell = 50; cell < 100; cell++)
  {
    waveProp.setHeight( cell,
                          0,
                          8899.326826472694 );
    waveProp.setMomentumX( cell,
                           0,
                           -122.0337839252433 );
	 waveProp.setBathymetry( cell,
	 								 0,
									 0 );
  }

  // set outflow boundary condition
  tsunami_lab::Boundary boundary[2] = {
		tsunami_lab::OUTFLOW, 
		tsunami_lab::OUTFLOW };
  waveProp.setGhostOutflow( boundary );

  // perform a time step
  for (int i = 0; i < 50; i++) {
		waveProp.setGhostOutflow( boundary );
		waveProp.timeStep(0.001, tsunami_lab::FWAVE );
  }

  // test for h*
  REQUIRE( waveProp.getHeight()[49]   == Approx(8899.739847378269) );
  REQUIRE( waveProp.getHeight()[50]   == Approx(8899.739847378269) );
}

TEST_CASE("Test the 1d wave propagation FWave solver (Rare-Rare Problem", "[WaveProp1dFWaveRareRare]")
{
  /*
   * @brief test state from middle_states.csv (Rare-Rare Problem)
   *
   * h_l = 9976.904476606509
   * h_r = 9976.904476606509
   * hu_l = -906.6229611756387
   * hu_r = 906.6229611756387
   * h* = 9974.006714260977
   */

  // construct solver and setup a Rare-Rare problem
  tsunami_lab::patches::WavePropagation1d waveProp( 100 );

  for (std::size_t cell = 0; cell < 50; cell++)
  {
    waveProp.setHeight( cell,
                          0,
                          9976.904476606509 );
    waveProp.setMomentumX( cell,
                             0,
                             -906.6229611756387 );
	 waveProp.setBathymetry( cell,
	 								 0,
									 0 );
  }
  for (std::size_t cell = 50; cell < 100; cell++)
  {
    waveProp.setHeight( cell,
                          0,
                          9976.904476606509 );
    waveProp.setMomentumX( cell,
                             0,
                             906.6229611756387 );
	 waveProp.setBathymetry( cell,
	 								 0,
									 0 );
  }

  // set outflow boundary condition
  tsunami_lab::Boundary boundary[2] = {
		tsunami_lab::OUTFLOW, 
		tsunami_lab::OUTFLOW };
  waveProp.setGhostOutflow( boundary );

  // perform a time step
  for (int i = 0; i < 50; i++) {
		waveProp.setGhostOutflow( boundary );
		waveProp.timeStep(0.001, tsunami_lab::FWAVE );
  }

  // test for h*
  REQUIRE( waveProp.getHeight()[49] == Approx(9974.006714260977) );
  REQUIRE( waveProp.getHeight()[50] == Approx(9974.006714260977) );
}