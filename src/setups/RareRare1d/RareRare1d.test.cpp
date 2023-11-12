/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Tests the rare-rare wave setup.
 **/
#include <catch2/catch.hpp>
#include "RareRare1d.h"

TEST_CASE( "Test the one-dimensional rare-rare wave setup.", "[RareRare1d]" ) {
  tsunami_lab::setups::RareRare1d l_RareRare( 25,
                                              55,
                                               3 );

  // left side
  REQUIRE( l_RareRare.getHeight( 2, 0 ) == 25 );

  REQUIRE( l_RareRare.getMomentumX( 2, 0 ) == 0 );

  REQUIRE( l_RareRare.getMomentumY( 2, 0 ) == 0 );

  REQUIRE( l_RareRare.getHeight( 2, 5 ) == 25 );

  REQUIRE( l_RareRare.getMomentumX( 2, 5 ) == 0 );

  REQUIRE( l_RareRare.getMomentumY( 2, 2 ) == 0 );

  // right side
  REQUIRE( l_RareRare.getHeight( 4, 0 ) == 55 );

  REQUIRE( l_RareRare.getMomentumX( 4, 0 ) == 0 );

  REQUIRE( l_RareRare.getMomentumY( 4, 0 ) == 0 );

  REQUIRE( l_RareRare.getHeight( 4, 5 ) == 55 );

  REQUIRE( l_RareRare.getMomentumX( 4, 5 ) == 0 );

  REQUIRE( l_RareRare.getMomentumY( 4, 2 ) == 0 );  
}