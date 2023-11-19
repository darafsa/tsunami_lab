/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Tests the subcritical case setup.
 **/
#include <catch2/catch.hpp>
#include "Subcritical1d.h"

TEST_CASE( "Test the one-dimensional subcritical case setup.", "[Subcritical1d]" ) {
  tsunami_lab::setups::Subcritical1d l_Subcritical;

  // x = 4: x between [0, 25] but below (8, 12)
  REQUIRE( l_Subcritical.getHeight( 4, 0 ) == 2 );

  REQUIRE( l_Subcritical.getMomentumX( 4, 0 ) == 4.42f );

  REQUIRE( l_Subcritical.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_Subcritical.getBathymetry( 4, 0 ) == -2 );

  // x = 10 (max Froude number): x between [0, 25] and (8, 12) 
  REQUIRE( l_Subcritical.getHeight( 10, 0 ) == 1.8f );

  REQUIRE( l_Subcritical.getMomentumX( 10, 0 ) == 4.42f );

  REQUIRE( l_Subcritical.getMomentumY( 10, 0 ) == 0 );

  REQUIRE( l_Subcritical.getBathymetry( 10, 0 ) == -1.8 );

  // x = 28 (max Froude number): x over [0, 25] and (8, 12) 
  REQUIRE( l_Subcritical.getHeight( 28, 0 ) == 0 );

  REQUIRE( l_Subcritical.getMomentumX( 28, 0 ) == 0 );

  REQUIRE( l_Subcritical.getMomentumY( 28, 0 ) == 0 );

  REQUIRE( l_Subcritical.getBathymetry( 28, 0 ) == -2 );
}