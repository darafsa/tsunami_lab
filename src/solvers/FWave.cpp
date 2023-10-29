/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * F-Wave solver for the shallow water equations.
**/

#include "FWave.h"

#include <cmath>

using namespace tsunami_lab::solvers;

void FWave::computeEigenvalues(float in_stateLeft[2], float in_stateRight[2], float out_eigenvaluesRoe[2]) {
	float heightLeft = in_stateLeft[0];
	float heightRight = in_stateRight[0];
	float momentumLeft = in_stateLeft[1];
	float momentumRight = in_stateRight[1];

	float sqrtHeightLeft = sqrt(heightLeft);
	float sqrtHeightRight = sqrt(heightRight);
	float particleVelocityLeft = momentumLeft / heightLeft;
	float particleVelocityRight = momentumRight / heightRight;
	
	float heightRoe = 0.5 * (heightLeft + heightRight);
	float particleVelocityRoe = particleVelocityLeft * sqrtHeightLeft + particleVelocityRight * sqrtHeightRight;
	particleVelocityRoe /= sqrtHeightLeft + sqrtHeightRight;

	float sqrtGTimesHeight = FWave::const_gSqrt * sqrt(heightRoe);
	
	out_eigenvaluesRoe[0] = particleVelocityRoe - sqrtGTimesHeight;
	out_eigenvaluesRoe[1] = particleVelocityRoe + sqrtGTimesHeight;
}

void FWave::computeInvertedEigenmatrix(float in_eigenvalues[2], float out_invertedEigenmatrix[2][2]) {
	float invertedMatrixDeterminant = 1 / (in_eigenvalues[1] - in_eigenvalues[0]);

	out_invertedEigenmatrix[0][0] =  invertedMatrixDeterminant * in_eigenvalues[1];
	out_invertedEigenmatrix[0][1] = -invertedMatrixDeterminant;
	out_invertedEigenmatrix[1][0] = -invertedMatrixDeterminant * in_eigenvalues[0];
	out_invertedEigenmatrix[1][1] =  invertedMatrixDeterminant;
}

void FWave::computeEigencoefficients(float in_stateLeft[2], float in_stateRight[2], float in_invertedEigenmatrix[2][2], float out_eigencoefficients[2]) {
	float heightLeft = in_stateLeft[0];
	float heightRight = in_stateRight[0];
	float momentumLeft = in_stateLeft[1];
	float momentumRight = in_stateRight[1];

	float flux_jump[2] = {
		momentumRight - momentumLeft,
		(momentumRight*momentumRight + 0.5f*FWave::const_g*heightRight*heightRight) - (momentumLeft*momentumLeft + 0.5f*FWave::const_g*heightLeft*heightLeft)
	};

	out_eigencoefficients[0] = in_invertedEigenmatrix[0][0] * flux_jump[0] + in_invertedEigenmatrix[0][1] * flux_jump[1];
	out_eigencoefficients[1] =	in_invertedEigenmatrix[1][0] * flux_jump[0] + in_invertedEigenmatrix[1][1] * flux_jump[1];
}

void FWave::computeNetUpdates(float in_stateLeft[2], float in_stateRight[2], float out_netUpdateLeft[2], float out_netUpdateRight[2]) {
	
	float eigenvalues[2];
	computeEigenvalues(in_stateLeft, in_stateRight, eigenvalues);

	float invertedEigenmatrix[2][2];
	computeInvertedEigenmatrix(eigenvalues, invertedEigenmatrix);

	float eigencoefficients[2];
	computeEigencoefficients(in_stateLeft, in_stateRight, invertedEigenmatrix, eigencoefficients);


	float waves[2][2] = {
		{ eigencoefficients[0], eigencoefficients[0] * eigenvalues[0] },
		{ eigencoefficients[1], eigencoefficients[1] * eigenvalues[1] }
	};

	out_netUpdateLeft[0] = 0;
	out_netUpdateLeft[1] = 0;
	out_netUpdateRight[0] = 0;
	out_netUpdateRight[1] = 0;

	for (int i = 0; i < 2; i++) {

		if( eigenvalues[i] < 0 ) {
			out_netUpdateLeft[0] += waves[i][0];
			out_netUpdateLeft[1] += waves[i][1];
		}
		else {
			out_netUpdateRight[0] += waves[i][0];
			out_netUpdateRight[1] += waves[i][1];
		}
  }
}
