/**
 * @author Marek Sommerfeld (marek.sommerfeld AT uni-jena.de)
 * @author Moritz Rätz (moritz.raetz AT uni-jena.de)
 *
 * @section DESCRIPTION
 * F-Wave solver for the shallow water equations.
**/

namespace tsunami_lab {
	namespace solvers {
		class FWave;
	}
}

class tsunami_lab::solvers::FWave {
	private:
		//! gravity constant
		static float constexpr const_g = 9.80665;
		//! squareroot of gravity constant
		static float constexpr const_gSqrt = 3.131557121;

		/**
		 * @brief Computes the Roe eigenvalues. 
		 * 
		 * @param in_stateLeft state of the left side; 0: height, 1: momentum.
		 * @param in_stateRight state of the right side; 0: height, 1: momentum.
		 * @param out_eigenvaluesRoe will be set to the Roe eigenvalues; 0: lambda^Roe_1, 1: lambda^Roe_2.
		 */
		static void computeEigenvalues(	float in_stateLeft[2], 
													float in_stateRight[2], 
													float out_eigenvaluesRoe[2] );

		/**
		 * @brief Calculates the flux function as a vaector.
		 * 
		 * @param in_state state of one cell; 0: height, 1: momentum.
		 * @param out_flux will be set to the flux function values; 0: hu, 1: hu^2+0.5gh^2.
		 */
		static void flux(float in_state[2], float out_flux[2]);

		/**
		 * @brief Computes the inverted eigenmatrix.
		 * 
		 * @param in_eigenvalues Roe eigenvalues; 0: lambda^Roe_1, 1: lambda^Roe_2.
		 * @param out_invertedEigenmatrix will be set to the inverted eigenmatrix; 0: eigenvector r1, 1: eigenvector r2.
		 */
		static void computeInvertedEigenmatrix(	float in_eigenvalues[2], 
																float out_invertedEigenmatrix[2][2] );

		/**
		 * @brief Computes the eigencoefficients.
		 * 
		 * @param in_stateLeft state of the left side; 0: height, 1: momentum.
		 * @param in_stateRight state of the right side; 0: height, 1: momentum.
		 * @param in_invertedEigenmatrix inverted eigenmatrix; 0: eigenvector r1, 1: eigenvector r2.
		 * @param out_eigencoefficients will be set to the eigencoefficients; 0: alpha_1, 1: alpha_2.
		 */
		static void computeEigencoefficients(	float in_stateLeft[2], 
															float in_stateRight[2], 
															float in_invertedEigenmatrix[2][2], 
															float out_eigencoefficients[2] );

	public:
		/**
		 * @brief Computes the net-updates.
		 * 
		 * @param in_stateLeft state of the left side; 0: height, 1: momentum.
		 * @param in_stateRight state of the right side; 0: height, 1: momentum.
		 * @param out_netUpdateLeft will be set to the net-updates for the left side; 0: height, 1: momentum.
		 * @param out_netUpdateRight will be set to the net-updates for the left side; 0: height, 1: momentum.
		 */
		static void computeNetUpdates(	float in_stateLeft[2], 
	 												float in_stateRight[2], 
													float out_netUpdateLeft[2], 
													float out_netUpdateRight[2] );
};
