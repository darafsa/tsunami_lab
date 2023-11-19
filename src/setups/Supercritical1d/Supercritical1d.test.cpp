/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Tests the supercritical case setup.
 **/
#include <catch2/catch.hpp>
#include "Supercritical1d.h"

TEST_CASE( "Test the one-dimensional supercritical case setup.", "[Supercritical1d]" ) {
  tsunami_lab::setups::Supercritical1d l_Supercritical;

  // x = 4: x between [0, 25] but below (8, 12)
  REQUIRE( l_Supercritical.getHeight( 4, 0 ) == 0.33f );

  REQUIRE( l_Supercritical.getMomentumX( 4, 0 ) == 0.18f );

  REQUIRE( l_Supercritical.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_Supercritical.getBathymetry( 4, 0 ) == -0.33f );

  // x = 10 (max Froude number): x between [0, 25] and (8, 12) 
  REQUIRE( l_Supercritical.getHeight( 10, 0 ) == 0.13f );

  REQUIRE( l_Supercritical.getMomentumX( 10, 0 ) == 0.18f );

  REQUIRE( l_Supercritical.getMomentumY( 10, 0 ) == 0 );

  REQUIRE( l_Supercritical.getBathymetry( 10, 0 ) == -0.13f );

  // x = 28 (max Froude number): x over [0, 25] and (8, 12) 
  REQUIRE( l_Supercritical.getHeight( 28, 0 ) == 0 );

  REQUIRE( l_Supercritical.getMomentumX( 28, 0 ) == 0 );

  REQUIRE( l_Supercritical.getMomentumY( 28, 0 ) == 0 );

  REQUIRE( l_Supercritical.getBathymetry( 28, 0 ) == -0.33f );
}