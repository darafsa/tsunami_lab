Riemann Solver | Project Report
===========================================================

Riemann Solver | Roe vs FWave
-------------------------------

**The FWave Solver is similar to the Roe Solver but differs in the following ways:**

instead of a jump in quantities

.. code-block:: c++

   // compute jump in quantities

  t_real l_hJump  = i_hR  - i_hL;
  t_real l_huJump = i_huR - i_huL;

a jump in fluxes is computed.

.. code-block:: c++

   float flux_jump[2] = {
		momentumRight - momentumLeft,
		(momentumRight*momentumRight + 0.5f*FWave::const_g*heightRight*heightRight) 
               - (momentumLeft*momentumLeft + 0.5f*FWave::const_g*heightLeft*heightLeft)
	};

instead of getting multiplied by :code:`l_sL` and :code:`l_sR`

.. code-block:: c++

   // compute scaled waves
  t_real l_waveL[2] = {0};
  t_real l_waveR[2] = {0};

  l_waveL[0] = l_sL * l_aL;
  l_waveL[1] = l_sL * l_aL * l_sL;

  l_waveR[0] = l_sR * l_aR;
  l_waveR[1] = l_sR * l_aR * l_sR;

the scaled waves are calculated without l_sL and l_sR (here eigencoefficients[0, 1]).

.. code-block:: c++

    float waves[2][2] = {
		{ eigencoefficients[0], eigencoefficients[0] * eigenvalues[0] },
		{ eigencoefficients[1], eigencoefficients[1] * eigenvalues[1] }
	};

instead of overriding each other

.. code-block:: c++

       // 1st wave
    if( l_sL < 0 ) {
      o_netUpdateL[l_qt] = l_waveL[l_qt];
    }
    else {
      o_netUpdateR[l_qt] = l_waveL[l_qt];
    }

    // 2nd wave
    if( l_sR > 0 ) {
      o_netUpdateR[l_qt] = l_waveR[l_qt];
    }
    else {
      o_netUpdateL[l_qt] = l_waveR[l_qt];

the netUpdates get added together (important regarding the supersonic problems).

.. code-block:: c++

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

Riemann Solver | Test Cases
---------------------------

Eigenvalues (Wave speeds):

.. code-block:: c++

  float stateLeft[2] =  {10, -3}
  float stateRight[2] = {9, 3};
  float eigenvaluesRoe[2]
  tsunami_lab::solvers::FWave::computeEigenvalues(stateLeft,
                                                stateLeft,
                                                eigenvaluesRoe);

  REQUIRE(eigenvaluesRoe[0] == Approx(-9.7311093998375095));
  EQUIRE(eigenvaluesRoe[1] == Approx(9.5731051658991654));

Inverted Eigenmatrix:

.. code-block:: c++

flux:

.. code-block:: c++

Eigencoefficients (Wave strengths):

.. code-block:: c++


Individual Member Contributions
--------------------------------

Marek Sommerfeld: Wrote base FWave Solver Code

Moritz Rätz: Wrote Unit test and Project Report


