2. Finite Volume Discretization | Project Report
================================================

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

   Bitte testcase einfügen

Implementation of shock-shock and rare-rare Problems
----------------------------------------------------

Implemented the shock-shock and rare-rare problems as setups. They are similar to the dam break structure, but differ in the water height and the momenta. 
There is only one water height, and the momentum is opposite for the two setups.

.. code-block:: c++

  Bitte Setup shock-shock und rare-rare einfügen

.. code-block:: c++

   Bitte Setup shock-shock und rare-rare einfügen

Observation of the influence of the initial parameters on shock-shock | rare-rare problems
------------------------------------------------------------------------------------------

Diskussion des Einflusses der Start Paramter :math:`h_l` und  :math:`u_l` und die connection von :math:`\lambda_{1/2} = u \mp \sqrt{gh}`
higher momenta -> higher middlestate height / does not affect wavespeed / higher initial height -> higher middlestate and faster wavespeeds 

something like this

.. math::
  h_r &= h_l \\
  hu_r &= -hu_l \\
  u_r &= -u_l \\
  h &= \frac{1}{2}(h_l+h_r) = h_l = h_r \\
  u &= \frac{u_l \sqrt{h_l} + u_r \sqrt{h_r}}{\sqrt{h_l}+\sqrt{h_r}} = \frac{u_l \sqrt{h} - u_l \sqrt{h}}{2\sqrt{h}} = 0 \\
  \lambda_{1,2} &= \mp \sqrt{gh}

Observation of the influence of the initial parameters on dam break problems
----------------------------------------------------------------------------

Diskussion über den Einfluss von initial water heights :math:`h_l` / :math:`h_r` und particle velocity :math:`u_r` auf fluss
higher initial height diff -> larger momentum in middle state / wavespeed of shock wave unaffected by right height  
higher inital height diff -> slower rare wave / wavespeed of shock shock proportional to sqrt of left height / higher momentum right side -> faster shockwave wavespeed (only minor)


Village Evacuation 
------------------

bla bla

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

Marek Sommerfeld: implementation of Solver interchangeability

Moritz Rätz: Wrote Unit test and Project Report


