/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Tests the shock-shock wave setup.
 **/
#include <catch2/catch.hpp>
#include "ShockShock1d.h"

TEST_CASE( "Test the one-dimensional shock-shock wave setup.", "[ShockShock1d]" ) {
  tsunami_lab::setups::ShockShock1d l_ShockShock( 25,
                                                  55,
                                                  3 );

  // left side
  REQUIRE( l_ShockShock.getHeight( 2, 0 ) == 25 );

  REQUIRE( l_ShockShock.getMomentumX( 2, 0 ) == 0 );

  REQUIRE( l_ShockShock.getMomentumY( 2, 0 ) == 0 );

  REQUIRE( l_ShockShock.getHeight( 2, 5 ) == 25 );

  REQUIRE( l_ShockShock.getMomentumX( 2, 5 ) == 0 );

  REQUIRE( l_ShockShock.getMomentumY( 2, 2 ) == 0 );

  // right side
  REQUIRE( l_ShockShock.getHeight( 4, 0 ) == 55 );

  REQUIRE( l_ShockShock.getMomentumX( 4, 0 ) == 0 );

  REQUIRE( l_ShockShock.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_ShockShock.getHeight( 4, 5 ) == 55 );

  REQUIRE( l_ShockShock.getMomentumX( 4, 5 ) == 0 );

  REQUIRE( l_ShockShock.getMomentumY( 4, 2 ) == 0 );  
}