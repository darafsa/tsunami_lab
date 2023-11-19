/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Tests the dam break setup.
 **/
#include <catch2/catch.hpp>
#include "DamBreak2d.h"

TEST_CASE( "Test the two-dimensional dam break setup.", "[DamBreak2d]" ) {
  tsunami_lab::setups::DamBreak2d l_damBreak( 5, 10, 10, 50, 50);

  // inner dam
  REQUIRE( l_damBreak.getHeight( 8, 4 ) == 10 );

  REQUIRE( l_damBreak.getMomentumX( 8, 4 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 8, 4 ) == 0 );

  REQUIRE( l_damBreak.getHeight( 3, 7 ) == 10 );

  REQUIRE( l_damBreak.getMomentumX( 3, 7 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 3, 7 ) == 0 );

  // outside of dam
  REQUIRE( l_damBreak.getHeight( 9, 7 ) == 5 );

  REQUIRE( l_damBreak.getMomentumX( 9, 7 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 9, 7 ) == 0 );

  REQUIRE( l_damBreak.getHeight( 12, 18 ) == 5 );

  REQUIRE( l_damBreak.getMomentumX( 12, 18 ) == 0 );

  REQUIRE( l_damBreak.getMomentumY( 12, 18 ) == 0 );
}