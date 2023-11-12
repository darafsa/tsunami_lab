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

void FWave::computeEigenvalues(real in_stateLeft[3], real in_stateRight[3], real out_eigenvaluesRoe[2]) {
	real heightLeft = in_stateLeft[0];
	real heightRight = in_stateRight[0];
	real momentumLeft = in_stateLeft[1];
	real momentumRight = in_stateRight[1];

	real sqrtHeightLeft = sqrt(heightLeft);
	real sqrtHeightRight = sqrt(heightRight);
	real particleVelocityLeft = momentumLeft / heightLeft;
	real particleVelocityRight = momentumRight / heightRight;
	
	real heightRoe = real(0.5) * (heightLeft + heightRight);
	real particleVelocityRoe = particleVelocityLeft * sqrtHeightLeft + particleVelocityRight * sqrtHeightRight;
	particleVelocityRoe /= sqrtHeightLeft + sqrtHeightRight;

	real sqrtGTimesHeight = const_gSqrt * sqrt(heightRoe);
	
	out_eigenvaluesRoe[0] = particleVelocityRoe - sqrtGTimesHeight;
	out_eigenvaluesRoe[1] = particleVelocityRoe + sqrtGTimesHeight;
}

void FWave::computeInvertedEigenmatrix(real in_eigenvalues[2], real out_invertedEigenmatrix[2][2]) {
	real invertedMatrixDeterminant = 1 / (in_eigenvalues[1] - in_eigenvalues[0]);

	out_invertedEigenmatrix[0][0] =  invertedMatrixDeterminant * in_eigenvalues[1];
	out_invertedEigenmatrix[0][1] = -invertedMatrixDeterminant;
	out_invertedEigenmatrix[1][0] = -invertedMatrixDeterminant * in_eigenvalues[0];
	out_invertedEigenmatrix[1][1] =  invertedMatrixDeterminant;
}

void FWave::flux(real in_state[3], real out_flux[2]) {
	real height = in_state[0];
	real momentum = in_state[1];

	out_flux[0] = momentum;
	out_flux[1] = (momentum * momentum / height + real(0.5) * const_g * height * height);
}

void FWave::computeDxPsi(real in_stateLeft[3], real in_stateRight[3], real & out_dxPsi) {
	real heightLeft = in_stateLeft[0];
	real heightRight = in_stateRight[0];
	real bathymetryLeft = in_stateLeft[2];
	real bathymetryRight = in_stateRight[2];
	
	//			    -g	  *	 bathymetryRight - bathymetryLeft  *  heightLeft + heightRight   / 2
	out_dxPsi = (-const_g * (bathymetryRight - bathymetryLeft) * ((heightLeft + heightRight) / 2));
}

void FWave::computeEigencoefficients(real in_stateLeft[3], real in_stateRight[3], real in_invertedEigenmatrix[2][2], real out_eigencoefficients[2]) {
	real fluxJumpLeft[2];
	real fluxJumpRight[2];
	real dxPsi;

	flux(in_stateLeft, fluxJumpLeft);
	flux(in_stateRight, fluxJumpRight);
	computeDxPsi(in_stateLeft, in_stateRight, dxPsi);
	
	real fluxJump[2] = {
		fluxJumpRight[0] - fluxJumpLeft[0],
		fluxJumpRight[1] - fluxJumpLeft[1] - dxPsi
	};

	out_eigencoefficients[0] = in_invertedEigenmatrix[0][0] * fluxJump[0] + in_invertedEigenmatrix[0][1] * fluxJump[1];
	out_eigencoefficients[1] =	in_invertedEigenmatrix[1][0] * fluxJump[0] + in_invertedEigenmatrix[1][1] * fluxJump[1];
}

void FWave::netUpdates(real in_stateLeft[3], real in_stateRight[3], real out_netUpdateLeft[2], real out_netUpdateRight[2]) {
	real eigenvalues[2];
	computeEigenvalues(in_stateLeft, in_stateRight, eigenvalues);

	real invertedEigenmatrix[2][2];
	computeInvertedEigenmatrix(eigenvalues, invertedEigenmatrix);

	real eigencoefficients[2];
	computeEigencoefficients(in_stateLeft, in_stateRight, invertedEigenmatrix, eigencoefficients);

	real waves[2][2] = {
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
