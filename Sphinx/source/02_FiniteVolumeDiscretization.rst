2. Finite Volume Discretization | Project Report
===========================================================

Integration of Interchangeability of Solver Type
------------------------------------------------

Altered the main() function to accept :code:`ROE` and :code:`FWAVE` as arguments.

.. code-block:: c++

  int main( int   i_argc,
            char *i_argv[] ) {
    // number of cells in x- and y-direction
    tsunami_lab::t_idx l_nx = 0;
    tsunami_lab::t_idx l_ny = 1;

    // set cell size
    tsunami_lab::t_real l_dxy = 1;

    // solver type
    tsunami_lab::patches::WavePropagation1d::Solver l_solverType;

    std::cout << "####################################" << std::endl;
    std::cout << "### Tsunami Lab                  ###" << std::endl;
    std::cout << "###                              ###" << std::endl;
    std::cout << "### https://scalable.uni-jena.de ###" << std::endl;
    std::cout << "####################################" << std::endl;

    if( i_argc != 3 ) {
      std::cerr << "invalid number of arguments, usage:" << std::endl;
      std::cerr << "  ./build/tsunami_lab N_CELLS_X SOLVER" << std::endl;
      std::cerr << "where N_CELLS_X is the number of cells in x-direction and SOLVER is the solver type to use." << std::endl;
      return EXIT_FAILURE;
    }
    else {
      l_nx = atoi( i_argv[1] );
      if( l_nx < 1 ) {
        std::cerr << "invalid number of cells" << std::endl;
        return EXIT_FAILURE;
      }
      l_dxy = 10.0 / l_nx;

    std::string l_solver = i_argv[2];
    if (l_solver == "FWAVE") {
      l_solverType = tsunami_lab::patches::WavePropagation::FWave;
    } else if (l_solver == "ROE") {
      l_solverType = tsunami_lab::patches::WavePropagation::Roe;
    } else {
      std::cerr << "invalid solver type. Please use either ROE or FWAVE" << std::endl;
        return EXIT_FAILURE;
    }
  }
    
    //(...)

    // iterate over time
    while( l_simTime < l_endTime ){
      
      //(...)
      
      l_waveProp->setGhostOutflow();
      std::string solver = "FWave";
      l_waveProp->timeStep( l_scaling, l_solverType );

      l_timeStep++;
      l_simTime += l_dt;
    }


Furthermore changed Wavepropagation1d.cpp for the same reason.

.. code-block:: c++

  #include "WavePropagation1d.h"
  #include "../solvers/FWave.h"
  #include "../solvers/Roe.h"
  
  void tsunami_lab::patches::WavePropagation1d::timeStep( t_real i_scaling, Solver i_solver ) {
    
    // pointers to old and new data
    t_real * l_hOld = m_h[m_step];
    t_real * l_huOld = m_hu[m_step];

    m_step = (m_step+1) % 2;
    t_real * l_hNew =  m_h[m_step];
    t_real * l_huNew = m_hu[m_step];

    // init new cell quantities
    for( t_idx l_ce = 1; l_ce < m_nCells+1; l_ce++ ) {
      l_hNew[l_ce] = l_hOld[l_ce];
      l_huNew[l_ce] = l_huOld[l_ce];
    }

    // iterate over edges and update with Riemann solutions
    for( t_idx l_ed = 0; l_ed < m_nCells+1; l_ed++ ) {
      // determine left and right cell-id
      t_idx l_ceL = l_ed;
      t_idx l_ceR = l_ed+1;

      // compute net-updates
      t_real l_netUpdates[2][2];
    
    t_real l_stateLeft[2] = {l_hOld[l_ceL], l_huOld[l_ceL]};
    t_real l_stateRight[2] = {l_hOld[l_ceR], l_huOld[l_ceR]};

      if ( i_solver == FWave ) {
        solvers::FWave::netUpdates( l_stateLeft, 
                          l_stateRight, 
                                    l_netUpdates[0],
                                    l_netUpdates[1] );
    } else {
      solvers::Roe::netUpdates( l_stateLeft[0], 
                        l_stateRight[0], 
                        l_stateLeft[1], 
                        l_stateRight[1], 
                                  l_netUpdates[0],
                                  l_netUpdates[1] );
      }

      // update the cells' quantities
      l_hNew[l_ceL]  -= i_scaling * l_netUpdates[0][0];
      l_huNew[l_ceL] -= i_scaling * l_netUpdates[0][1];

      l_hNew[l_ceR]  -= i_scaling * l_netUpdates[1][0];
      l_huNew[l_ceR] -= i_scaling * l_netUpdates[1][1];
    }
  }

Usage of the middle states csv as sanity check
----------------------------------------------

Added new test cases in the Wavepropagation1d unit tests for the rare-rare and shock-shock problems.

(Example for test case of rare-rare problem)

.. code-block:: c++

  TEST_CASE("Test the 1d wave propagation FWave solver (Rare-Rare Problem", "[WaveProp1dFWaveRareRare]")
  {
  /**
   * @brief test state from middle_states.csv (Rare-Rare Problem)
   *
   * h_l = 7589.71304876485
   * h_r = 7589.71304876485
   * hu_l = -138.9853242339589
   * hu_r = 138.9853242339589
   * h* = 7589.203700916305
   */

  // construct solver and setup a Rare-Rare problem
  tsunami_lab::patches::WavePropagation1d m_waveProp(100);

  for (std::size_t l_ce = 0; l_ce < 50; l_ce++)
    {
    m_waveProp.setHeight(l_ce,
                         0,
                         7589.71304876485);
    m_waveProp.setMomentumX(l_ce,
                            0,
                            -138.9853242339589);
    }
  for (std::size_t l_ce = 50; l_ce < 100; l_ce++)
    {
    m_waveProp.setHeight(l_ce,
                         0,
                         7589.71304876485);
    m_waveProp.setMomentumX(l_ce,
                            0,
                            138.9853242339589);
    }

    // set outflow boundary condition
    m_waveProp.setGhostOutflow();

    // perform a time step
    for (int i = 0; i < 30; i++)
      {
        m_waveProp.timeStep(0.001);
      }

    // test for h*
    REQUIRE(m_waveProp.getHeight()[49] == Approx(7589.203700916305));
    REQUIRE(m_waveProp.getHeight()[50] == Approx(7589.203700916305));
  }

Implementation of shock-shock and rare-rare Problems
----------------------------------------------------

Implemented the shock-shock and rare-rare problems as setups. They are similar to the dam break structure, but differ in the water height and the momenta. 
There is only one water height, and the momentum is opposite for the two setups.

.. code-block:: c++

  tsunami_lab::setups::ShockShock1d::ShockShock1d(t_real i_height,
                                                  t_real i_momentum,
                                                  t_real i_midPos)
  {
      m_height = i_height;
      m_momentum = i_momentum;
      m_midPos = i_midPos;
  }

 tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getHeight(t_real,
                                                                  t_real) const
  {
    return m_height;
  }

  tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getMomentumX(t_real i_x,
                                                                      t_real) const
  {
    if (i_x < m_midPos)
      {
          return m_momentum;
      }
    else
      {
          return -m_momentum;
      }
  }

  tsunami_lab::t_real tsunami_lab::setups::ShockShock1d::getMomentumY(t_real,
                                                                      t_real) const
  {
    return 0;
  }

.. code-block:: c++

  tsunami_lab::setups::RareRare1d::RareRare1d(t_real i_height,
                                              t_real i_momentum,
                                              t_real i_midPos)
  {
      m_height = i_height;
      m_momentum = i_momentum;
      m_midPos = i_midPos;
  }

  tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getHeight(t_real,
                                                                t_real) const
  {
    return m_height;
  }

  tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getMomentumX(t_real i_x,
                                                                    t_real) const
  {
    if (i_x < m_midPos)
      {
          return -m_momentum;
      }
    else
      {
          return m_momentum;
      }
  }

  tsunami_lab::t_real tsunami_lab::setups::RareRare1d::getMomentumY(t_real,
                                                                    t_real) const
  {
    return 0;
  }

Observation of the influence of initial parameters on shock-shock/rare-rare problems
------------------------------------------------------------------------------------

The higher the momentum :math:`u_l`, the higher the middlestate height.
The higher the height :math:`h_l`, the higher the middlestate height and wavespeeds.

.. video:: _static/animations/02/rare_10_5_70.mp4
	 :autoplay:
	 :nocontrols:
	 :loop:
	 :height: 300
	 :width: 650

.. video:: _static/animations/02/rare_10_10_70.mp4
	 :autoplay:
	 :nocontrols:
	 :loop:
	 :height: 300
	 :width: 650

The Wavespeeds are proportional to the square root of the initial height. Furthermore they are independent of the initial momentum.

.. math::
  h_r &= h_l \\
  hu_r &= -hu_l \\
  u_r &= -u_l \\
  h &= \frac{1}{2}(h_l+h_r) = h_l = h_r \\
  u &= \frac{u_l \sqrt{h_l} + u_r \sqrt{h_r}}{\sqrt{h_l}+\sqrt{h_r}} = 0 \\
  \lambda_{1,2} &= \mp \sqrt{gh}

Observation of the influence of the initial parameters on dam break problems
----------------------------------------------------------------------------

The higher the initial height difference, the higher the middle state. It affects the rarefaction wave as well in terms of wave speed (slowing down). 
Furthermore, it seems that the (shock) wavespeeds are proportional to the sqrt of the height.


Village Evacuation 
------------------

initial values:

.. math::

  s_{village} = 25km = 25000m\quad q_l = \begin{bmatrix} 14 \\ 0 \end{bmatrix}\quad q_r = \begin{bmatrix} 3.5 \\ 0.7 \end{bmatrix}\\

with

.. math::

  h_l = 14m\quad h_r = 3.5m\quad u_l = 0 \frac{m}{s}\quad u_r = 0.7 \frac{m}{s}

calculate Roe height :math:`h^{Roe}` and the Roe wavespeed :math:`u^{Roe}`:

.. math::

  h^{Roe} &= \frac{1}{2} (h_l + h_r) = \frac{1}{2} (14m + 3.5m) = 8.75 m \\
  u^{Roe} &= \frac{u_l \sqrt{h_l} + u_r \sqrt{h_r}}{\sqrt{h_l}+\sqrt{h_r}} = \frac{0 \frac{m}{s} \cdot \sqrt{14m} + 0.7 \frac{m}{s} \cdot \sqrt{3.5m}}{\sqrt{14m}+\sqrt{3.5m}} = 0.23333 \frac{m}{s}\\

now use :math:`h^{Roe}` and :math:`u^{Roe}` to calculate :math:`\lambda_r^{Roe}`:

.. math::

  \lambda_r^{Roe} = u^{Roe} + \sqrt{gh^{Roe}} = 0.23333\frac{m}{s} + \sqrt{9.80665\frac{m}{s^2} \cdot 8.75m} = 9.49660 \frac{m}{s} \\

and finally calculate the time left for evacuation of the village :math:`t_{evacuation}`:

.. math::

  t_{evacuation} = \frac{s_{village}}{\lambda_r^{Roe}} = \frac{25000m}{9.49660 \frac{m}{s}} = 2.632,52 s = 43.88 min


Individual Member Contributions
--------------------------------

Marek Sommerfeld: implementation of Solver interchangeability, Rare-Rare/Shock-Shock setup, CI Integration unit tests/Sphinx static page

Moritz Rätz: Wavepropagation Unit test, Rare-Rare/Shock-Shock setup and Project Report


