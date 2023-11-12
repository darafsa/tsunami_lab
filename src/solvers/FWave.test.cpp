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
 *
 *  h: 10 | 9
 * hu: -3 | 3
 *  b:  0 | 0
 *
 * FWave height: 9.5
 * FWave velocity: ((-3/10) * sqrt(10) + (3/9) * sqrt(9)) / (sqrt(10) + sqrt(9))
 *                  = -0.0083275543199207307978
 * FWave Eigenvalues (speeds): s1 = velocity - sqrt(9.80665 * height) = -9.64378
 *                             s2 = velocity + sqrt(9.80665 * height) =  9.66043
 *
 * wolframalpha.com query: ((-3/10) * sqrt(10) + (3/9) * sqrt(9)) / (sqrt(10) + sqrt(9)) - sqrt(9.80665 * 9.5)
 */

float stateLeft[3] =  {10, -3, 0};
float stateRight[3] = {9, 3, 0};
float eigenvaluesRoe[2];
tsunami_lab::solvers::FWave::computeEigenvalues( stateLeft,
                                                 stateRight,
                                                 eigenvaluesRoe );

REQUIRE(eigenvaluesRoe[0] == Approx(-9.64378));
REQUIRE(eigenvaluesRoe[1] == Approx(9.66043));
}

TEST_CASE("Test the computation of the InvertedEigenmatrix.", "[FWaveInvertedEigenmatrix]")
{
/*
 * Test case:
 *
 *  h:  10 | 9
 *  u:  -3 | 3
 * hu: -30 | 27
 *  b:   0 | 0
 *
 * The derivation of the Eigenvalues is given above.
 *
 *  Matrix of Eigenvectors:
 *
 *      | 1   1 |
 *  R = |       |
 *      | s1 s2 |
 *
 * Inversion yields:
 * 
 * wolframalpha.com query: invert {{1, 1}, {-9.7311093998375095, 9.5731051658991654}}
 *
 *        | 0.49590751974393229 -0.051802159398648326 |
 * Rinv = |                                           |
 *        | 0.50409248025606771  0.051802159398648326 |
 *
 */

float eigenvalues[2] = {-9.7311093998375095, 9.5731051658991654};
float invertedEigenmatrix[2][2];
tsunami_lab::solvers::FWave::computeInvertedEigenmatrix( eigenvalues,
                                                         invertedEigenmatrix );

REQUIRE(invertedEigenmatrix[0][0] == Approx(0.49590751974393229));
REQUIRE(invertedEigenmatrix[0][1] == Approx(-0.051802159398648326));
REQUIRE(invertedEigenmatrix[1][0] == Approx(0.50409248025606771));
REQUIRE(invertedEigenmatrix[1][1] == Approx(0.051802159398648326));
}

TEST_CASE("Test the computation of the Eigencoefficients (wave strenghts).", "[FWaveEigencoefficients]")
{
/*
 * Test case:
 *
 *  h:  10 | 9
 *  u:  -3 | 3
 * hu: -30 | 27
 *  b:   0 | 0
 *
 * The derivation of the Eigenmatrix is given above.
 *
 * Multiplicaton with the jump in fluxes gives the Eigencoefficients:
 *
 * wolframalpha.com query: {{0.49590751974393229, -0.051802159398648326}, {0.50409248025606771, 0.051802159398648326}} * {27--30, 9*3^2+1/2*9.80665*9^2-(10*(-3)^2+1/2*9.80665*10^2)}
 *
 *  | 0.49590751974393229 -0.051802159398648326 |   |       27              -         -30                |   |   33.559  |
 *  |                                           | * |                                                    | = |           |
 *  | 0.50409248025606771  0.051802159398648326 |   | 9*3^2+1/2*9.80665*9^2 - (10*-3^2+1/2*9.80665*10^2) |   |   23.44   |
 *
 */

float eigenmatrix[2][2] = {{0.49590751974393229, -0.051802159398648326}, {0.50409248025606771, 0.051802159398648326}};
float stateLeft[3] = {10, -30, 0};
float stateRight[3] = {9, 27, 0};
float eigencoefficients[2];
tsunami_lab::solvers::FWave::computeEigencoefficients( stateLeft,
                                                       stateRight,
                                                       eigenmatrix,
                                                       eigencoefficients );

REQUIRE(eigencoefficients[0] == Approx(33.559));
REQUIRE(eigencoefficients[1] == Approx(23.441));
}


TEST_CASE("Test the derivation of the FWave net-updates.", "[FWaveUpdates]")
{
/*
 * Test case:
 *
 *      left | right
 *  h:    10 | 9
 *  u:    -3 | 3
 * hu:   -30 | 27
 *  b:     0 | 0
 *
 * The derivation of the FWave Eigenvalues (s1, s2) and eigencoefficients (a1, a1) is given above.
 *
 * The net-updates are given through the scaled eigenvectors.
 *
 *                      |  1 |   |         33.559            |
 * update #1:      a1 * |    | = |                           |
 *                      | s1 |   | -326.5663003491469813105  |
 *
 *                      |  1 |   |          23.441           |
 * update #2:      a2 * |    | = |                           |
 *                      | s2 |   | 224.4031581938423361414   |
 */

float stateLeft[3] = {10, -30, 0};
float stateRight[3] = {9, 27, 0};
float netUpdateLeft[2];
float netUpdateRight[2];

tsunami_lab::solvers::FWave::netUpdates( stateLeft,
                                         stateRight,
                                         netUpdateLeft,
                                         netUpdateRight );

REQUIRE(netUpdateLeft[0] == Approx(33.559));
REQUIRE(netUpdateLeft[1] == Approx(-326.5663003491469813105));

REQUIRE(netUpdateRight[0] == Approx(23.441));
REQUIRE(netUpdateRight[1] == Approx(224.4031581938423361414));

/*
* Test case (dam break):
*
*     left | right
*   h:  10 | 8
*  hu:   0 | 0
*   b:   0 | 0
*
* FWave speeds are given as:
*
*   s1 = -sqrt(9.80665 * 9)
*   s2 =  sqrt(9.80665 * 9)
*
* Inversion of the matrix of right Eigenvectors:
*
*   wolframalpha.com query: invert {{1, 1}, {-sqrt(9.80665 * 9), sqrt(9.80665 * 9)}}
*
*          | 0.5 -0.0532217 |
*   Rinv = |                |
*          | 0.5 -0.0532217 |
*
*   wolframalpha.com query: {{0.5,-0.0532217},{0.5,0.0532217}} * {{0},{1/2*9.80665*8^2 - 1/2*9.80665*10^2}}
*
* Multiplicaton with the jump in fluxes gives the wave strengths:
*
*        |  0              -          0       |   |  9.39468  |   | a1 |
* Rinv * |                                    | = |           | = |    |
*        | 1/2*9.80665*8^2 - 1/2*9.80665*10^2 |   | -9.39468  |   | a2 |
*
* The net-updates are given through the scaled eigenvectors.
*
*                      |  1 |   |    9.39468    |
* update #1:      a1 * |    | = |               |
*                      | s1 |   |   -88.2599    |
*
*                      |  1 |   |   -9.39468    |
* update #2:      a2 * |    | = |               |
*                      | s2 |   |   -88.2599    |
*/

stateLeft[0] = 10;
stateLeft[1] = 0;
stateRight[0] = 8;
stateRight[1] = 0;

tsunami_lab::solvers::FWave::netUpdates( stateLeft,
                                         stateRight,
                                         netUpdateLeft,
                                         netUpdateRight );

REQUIRE(netUpdateLeft[0] == Approx(9.39468));
REQUIRE(netUpdateLeft[1] == Approx(-88.2599));

REQUIRE(netUpdateRight[0] == Approx(-9.39468));
REQUIRE(netUpdateRight[1] == Approx(-88.2599));

/*
 * Test case supersonic problem
 *
 *      left | right
 *  h:    1  | 1
 *  u:   100 | 10
 * hu:   100 | 10
 *  b:    0  | 0
 *
 *  FWave Eigenvalues are given as:
 *
 *  s1 = 55 - sqrt(9.80665 * 1) = 51.868443
 *  s2 = 55 + sqrt(9.80665 * 1) = 58.131557
 *
 *  Inversion of the Eigenmatrix:
 *
 *   wolframalpha.com query: invert {{1, 1}, {55 - sqrt(9.80665 * 1), 55 + sqrt(9.80665 * 1)}}
 *
 *          | 9.28157  -0.159665  |
 *   Rinv = |                     |
 *          | -8.28157  0.159665  |
 *
 * Multiplicaton with the jump in fluxes gives the eigencoefficients:
 *
 *        |        10             -          100             |   |  745.342  |   | a1 |
 * Rinv * |                                                  | = |           | = |    |
 *        |1*10^2+1/2*9.80665*1^2 - (1*100^2+1/2*9.80665*1^2)|   | -835.342  |   | a2 |
 *
 *
 * update #1:     0
 *
 *                     |  1 |         |  1 |    |      -90         |
 * update #2:     a1 * |    | +  a2 * |    |  = |                  |
 *                     | s1 |         | s2 |    | -9900.002044988  |
 */

stateLeft[0] = 1;
stateLeft[1] = 100;
stateRight[0] = 1;
stateRight[1] = 10;

tsunami_lab::solvers::FWave::netUpdates( stateLeft,
                                         stateRight,
                                         netUpdateLeft,
                                         netUpdateRight );

REQUIRE(netUpdateLeft[0] == Approx(0));
REQUIRE(netUpdateLeft[1] == Approx(0));

REQUIRE(netUpdateRight[0] == Approx(-90));
REQUIRE(netUpdateRight[1] == Approx(-9900.002044988));

/*
* Test case (dam break) with bathymetry:
*
*     left | right
*   h:  10 |  8
*  hu:   0 |  0
*   b:  -6 | -13
*
* FWave speeds are given as:
*
*   s1 = -sqrt(9.80665 * 9)
*   s2 =  sqrt(9.80665 * 9)
*
* Inversion of the matrix of right Eigenvectors:
*
*   wolframalpha.com query: invert {{1, 1}, {-sqrt(9.80665 * 9), sqrt(9.80665 * 9)}}
*
*          | 0.5    -0.0532217 |
*   Rinv = |                   |
*          | 0.5    -0.0532217 |
*
* new bathymetry (dxPsi) computation:
*
* -g * (bathymetryRight - bathymetryLeft) * ((heightLeft + heightRight) / 2)
*
* wolframalpha.com query: {-9.80665 * (-13 + 6) * ((10 + 8) / 2)}
* 
* -9.80665 * (-13 + 6) * ((10 + 8) / 2) = 617.81895
*
* Multiplicaton with the jump in fluxes gives the wave strengths:
*
* wolframalpha.com query: {{0.5,-0.0532217},{0.5,0.0532217}} * {{0},{1/2*9.80665*8^2 - 1/2*9.80665*10^2 + 9.80665 * (-13 + 6) * (10 + 8) / 2}}
*
*        |  0              -          0                   |   |  42.2761  |   | a1 |
* Rinv * |                                                | = |           | = |    |
*        | 1/2*9.80665*8^2 - 1/2*9.80665*10^2 - 617.81895 |   | -42.2761  |   | a2 |
*
* The net-updates are given through the scaled eigenvectors.
*
* wolframalpha.com query: 42.2761* {{1},{-sqrt(9.80665 * 9)}} for update #1
*
*                      |  1 |   |   42.2761   |
* update #1:      a1 * |    | = |             |
*                      | s1 |   |  -397.17    |
*
* wolframalpha.com query: -42.2761* {{1},{sqrt(9.80665 * 9)}} for update #2
*
*                      |  1 |   |  -42.2761   |
* update #2:      a2 * |    | = |             |
*                      | s2 |   |  -397.17    |
*/

stateLeft[0] = 10;
stateLeft[1] = 0;
stateLeft[2] = -6;
stateRight[0] = 8;
stateRight[1] = 0;
stateRight[2] = -13;

tsunami_lab::solvers::FWave::netUpdates( stateLeft,
                                         stateRight,
                                         netUpdateLeft,
                                         netUpdateRight );

REQUIRE(netUpdateLeft[0] == Approx(42.2761));
REQUIRE(netUpdateLeft[1] == Approx(-397.17));

REQUIRE(netUpdateRight[0] == Approx(-42.2761));
REQUIRE(netUpdateRight[1] == Approx(-397.17));
}