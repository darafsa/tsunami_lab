/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Unit tests of the F-Wave solver.
 **/
#include <catch2/catch.hpp>
#define private public
#include "FWave.h"
#undef public

TEST_CASE("Test the computation of the Eigenvalues (FWave speeds).", "[FWaveEigenvalues]")
{
/*
 * Test case:
 *  h: 10 | 9
 *  u: -3 | 3
 *
 * FWave height: 9.5
 * FWave velocity: (sqrt(10) * -3 + 3 * 3) / ( sqrt(10) + sqrt(9) )
 *                  = -0.0790021169691720
 * FWave speeds: s1 = -0.079002116969172024 - sqrt(9.80665 * 9.5) = -9.7311093998375095
 *               s2 = -0.079002116969172024 + sqrt(9.80665 * 9.5) =  9.5731051658991654
 */

float stateLeft[2] =  {10, -3}
float stateRight[2] = {9, 3};
float eigenvaluesRoe[2]
tsunami_lab::solvers::FWave::computeEigenvalues(stateLeft,
                                                stateLeft,
                                                eigenvaluesRoe);

REQUIRE(eigenvaluesRoe[0] == Approx(-9.7311093998375095));
REQUIRE(eigenvaluesRoe[1] == Approx(9.5731051658991654));
}

TEST_CASE("Test the computation of the InvertedEigenmatrix.", "[FWaveInvertedEigenmatrix]")
{
/*
 * Test case:
 * Eigenvalues: 1 | 5
 *
 * invertedMatrixDeterminant: 1 / (1 - 5) = -1/4
 *
 *        | -1/4 * 1        -(-1/4) |
 * Rinv = |                         |
 *        | -(-1/4) * 1        -1/4 |
 *
 *        | -1/4                1/4 |
 * Rinv = |                         |
 *        | 1/4                -1/4 |
 */

float eigenvalues[2] = {1, 5}
float invertedEigenmatrix[2][2]
tsunami_lab::solvers::FWave::computeInvertedEigenmatrix(eigenvalues,
                                                        invertedEigenmatrix);

REQUIRE(invertedEigenmatrix[0][0] == Approx(-0.25));
REQUIRE(invertedEigenmatrix[0][1] == Approx(0.25));
REQUIRE(invertedEigenmatrix[1][0] == Approx(0.25));
REQUIRE(invertedEigenmatrix[1][1] == Approx(-0.25));
}


TEST_CASE("Test the computation of the FWave Eigencoefficients.", "[FWaveEigencoefficients]")
{
/*
 *
 */


REQUIRE( == Approx());
REQUIRE( == Approx());
}

TEST_CASE("Test the derivation of the FWave net-updates.", "[FWaveUpdates]")
{
/*
 * Test case
 *
 *      left | right
 *  h:    5  | 3
 *  u:   -6  | 9
 *  hu:  -15 | 21
 *
 */

float stateLeft[2] = {5, -15};
float stateRight[2] = {3, 21};
float netUpdateLeft[2];
float netUpdateRight[2];
tsunami_lab::solvers::FWave::netUpdates(stateLeft,
                                        stateRight,
                                        netUpdateLeft,
                                        netUpdateRight);

REQUIRE(netUpdateLeft[0] == Approx(10.942));
REQUIRE(netUpdateLeft[1] == Approx(-53.5962));

REQUIRE(netUpdateRight[0] == Approx(25.058));
REQUIRE(netUpdateRight[1] == Approx(191.143));

/*
 * Test case dam break
 *
 *      left | right
 *  h:   13  | 11
 *  hu:    0 | 0
 *
 */

float stateLeft[2]  = {13, 0};
float stateRight[2] = {11, 0};
float netUpdateLeft[2];
float netUpdateRight[2];
tsunami_lab::solvers::FWave::netUpdates(stateLeft,
                                        stateRight,
                                        netUpdateLeft,
                                        netUpdateRight);

REQUIRE(netUpdatesLeft[0] == Approx(?));
REQUIRE(netUpdatesLeft[1] == Approx(?));

REQUIRE(netUpdatesRight[0] == Approx(?));
REQUIRE(netUpdatesRight[1] == Approx(?));

/*
 * Test case supersonic problem
 *
 *      left | right
 *  h:   1   | 1
 *  u:   100 | 10
 *  hu:  100 | 100
 *
 */

float stateLeft[2] =  {1, 100};
float stateRight[2] = {1, 10};
float netUpdateLeft[2];
float netUpdateRight[2];
tsunami_lab::solvers::FWave::netUpdates(stateLeft,
                                        stateRight,
                                        netUpdateLeft,
                                        netUpdateRight);

REQUIRE(netUpdatesLeft[0] == Approx(0));
REQUIRE(netUpdatesLeft[1] == Approx(0));

REQUIRE(netUpdatesRight[0] == Approx(-90));
REQUIRE(netUpdatesRight[1] == Approx(-9900.002044988));
