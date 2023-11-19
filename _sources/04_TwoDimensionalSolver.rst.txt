4. Two Dimensional Solver | Project Report
===========================================================

Dimensional Splitting
---------------------

We achieve two-dimensionality through dimensional splitting, i.e. first the x-sweep is calculated where: 

:math:`Q_{i,j}^* = Q_{i,j}^n - \frac{\Delta t}{\Delta x} \left( A^+ \Delta Q_{i-1/2,j} + A^- \Delta Q_{i+1/2,j} \right)  \quad \forall i \in \{ 1, .., n \}, \; j \in \{ 0, .., n+1 \}`

and then the y-sweep is calculated by using the result of the x-sweep :math:`Q_{i,j}^*` :

:math:`Q_{i,j}^{n+1} = Q_{i,j}^* - \frac{\Delta t}{\Delta y} \left( B^+ \Delta Q^*_{i,j-1/2} + B^- \Delta Q^*_{i,j+1/2} \right)  \quad \forall i,j \in \{ 1, .., n \}`

This was implemented through a new Wavepropagation2d patch to preserve the functionality of the one dimensional solver. 

Now there needs to be more memory allocated and deleted: 

.. code-block:: c++
  
  WavePropagation2d::WavePropagation2d( idx in_cellCountX, idx in_cellCountY ) {
    cellCountX = in_cellCountX;
    cellCountY = in_cellCountY;

    // allocate memory including a single ghost cell on each side
    for( unsigned short step = 0; step < 2; step++ ) {
      height[step] = new real*[ cellCountX + 2 ];
      momentumX[step] = new real*[ cellCountX + 2 ];
      momentumY[step] = new real*[ cellCountX + 2 ];
      for ( idx x = 0; x < cellCountX + 2; x++) {
        height[step][x] = new real[ cellCountY + 2];
        momentumX[step][x] = new real[ cellCountY + 2];
        momentumY[step][x] = new real[ cellCountY + 2];
      }
    }
    bathymetry = new real*[ cellCountX + 2 ];
    for ( idx x = 0; x < cellCountX + 2; x++) {
      bathymetry[x] = new real[ cellCountY + 2];
    }
    array1d = new real[ cellCountX * cellCountY ];

    // init to zero
    for( unsigned short step = 0; step < 2; step++ ) {
      for( idx y = 0; y < cellCountY + 2; y++ ) {
        for( idx x = 0; x < cellCountX + 2; x++ ) {
          height[step][x][y] = 0;
          momentumX[step][x][y] = 0;
          momentumY[step][x][y] = 0;
          if( step == 0 ) {
            bathymetry[x][y] = 0;
          }
        }
      }
    }
    for ( idx cell = 0; cell < cellCountX * cellCountY; cell++ ) {
      array1d[cell] = 0;
    }
  }

  WavePropagation2d::~WavePropagation2d() {
    for( unsigned short step = 0; step < 2; step++ ) {
      for ( idx x = 0; x < cellCountX + 2; x++) {
        delete height[step][x];
        delete momentumX[step][x];
        delete momentumY[step][x];
      }
      delete[] height[step];
      delete[] momentumX[step];
      delete[] momentumY[step];
    }
    for ( idx x = 0; x < cellCountX + 2; x++) {
      delete bathymetry[x];
    }
    delete[] bathymetry;
  }

Timestep function x-sweep:

.. code-block:: c++

  void WavePropagation2d::timeStep( real in_scaling, Solver in_solver ) {
    // pointers to old and new data
    real ** heightOld = height[step];
    real ** momentumXOld = momentumX[step];
    real ** momentumYOld = momentumY[step];

    step = (step+1) % 2;
    real ** heightNew =	height[step];
    real ** momentumXNew = momentumX[step];
    real ** momentumYNew = momentumY[step];

    // init new cell quantities
    for( idx y = 1; y < cellCountY + 1; y++ ) {
      for( idx x = 1; x < cellCountX + 1; x++) {
        heightNew[x][y] = heightOld[x][y];
        momentumXNew[x][y] = momentumXOld[x][y];
        momentumYNew[x][y] = momentumYOld[x][y];
      }
    }

    // iterate over edges and update with Riemann solutions in x-direction
    for( idx y = 0; y < cellCountY + 2; y++ ) {
      for( idx edgeX = 0; edgeX < cellCountX + 1; edgeX++ ) {
        // determine cell-id
        idx cellLeft = edgeX;
        idx cellRight = edgeX+1;

        // compute net-updates
        real netUpdates[2][2];
      
        real stateLeft[3] = { heightOld[cellLeft][y], momentumXOld[cellLeft][y], bathymetry[cellLeft][y] };
        real stateRight[3] = { heightOld[cellRight][y], momentumXOld[cellRight][y], bathymetry[cellRight][y] };

        if(bathymetry[cellLeft][y] > 0) {
          stateLeft[0] = stateRight[0];
          stateLeft[1] = -stateRight[1];
          stateLeft[2] = stateRight[2];
        }

        if(bathymetry[cellRight][y] > 0) {
          stateRight[0] = stateLeft[0];
          stateRight[1] = -stateLeft[1];
          stateRight[2] = stateLeft[2];
        }

        if ( in_solver == FWAVE ) {
          solvers::FWave::netUpdates( stateLeft, stateRight, netUpdates[0], netUpdates[1] );
        } else {
          solvers::Roe::netUpdates( stateLeft[0], stateRight[0], stateLeft[1], stateRight[1], netUpdates[0], netUpdates[1] );
        }

        // update the cells' quantities
        heightNew[cellLeft][y] -= in_scaling * netUpdates[0][0];
        momentumXNew[cellLeft][y] -= in_scaling * netUpdates[0][1];

        heightNew[cellRight][y]	-= in_scaling * netUpdates[1][0];
        momentumXNew[cellRight][y] -= in_scaling * netUpdates[1][1];
      }
    }

and then the y-sweep:

.. code-block:: c++

  // iterate over edges and update with Riemann solutions in y-direction
    for( idx x = 0; x < cellCountX + 2; x++ ) {
      for( idx edgeY = 0; edgeY < cellCountX + 1; edgeY++ ) {
        // determine cell-id
        idx cellTop = edgeY+1;
        idx cellBottom = edgeY;

        // compute net-updates
        real netUpdates[2][2];
      
        real stateLeft[3] = { heightOld[x][cellBottom], momentumYOld[x][cellBottom], bathymetry[x][cellBottom] };
        real stateRight[3] = { heightOld[x][cellTop], momentumYOld[x][cellTop], bathymetry[x][cellTop] };

        if(bathymetry[x][cellBottom] > 0) {
          stateLeft[0] = stateRight[0];
          stateLeft[1] = -stateRight[1];
          stateLeft[2] = stateRight[2];
        }

        if(bathymetry[x][cellTop] > 0) {
          stateRight[0] = stateLeft[0];
          stateRight[1] = -stateLeft[1];
          stateRight[2] = stateLeft[2];
        }

        if ( in_solver == FWAVE ) {
          solvers::FWave::netUpdates( stateLeft, stateRight, netUpdates[0], netUpdates[1] );
        } else {
          solvers::Roe::netUpdates( stateLeft[0], stateRight[0], stateLeft[1], stateRight[1], netUpdates[0], netUpdates[1] );
        }

        // update the cells' quantities
        heightNew[x][cellBottom] -= in_scaling * netUpdates[0][0];
        momentumYNew[x][cellBottom] -= in_scaling * netUpdates[0][1];

        heightNew[x][cellTop] -= in_scaling * netUpdates[1][0];
        momentumYNew[x][cellTop] -= in_scaling * netUpdates[1][1];
      }
    }

For the ghost cells a copyGhostCells function was implemented (respectively for Outflow and Recflection). Therefore, the original setGhostOutflow function was altered:

.. code-block:: c++

  void WavePropagation2d::setGhostOutflow( Boundary in_boundary[2] ) {
    // set left boundary
    if(in_boundary[0] == OUTFLOW) {
      copyGhostCellsOutflow( height[step] );
      copyGhostCellsOutflow( momentumX[step] );
      copyGhostCellsOutflow( momentumY[step] );
      copyGhostCellsOutflow( bathymetry );
    } else if (in_boundary[0] == REFLECTING) {
      copyGhostCellsReflecting( height[step], 0 );
      copyGhostCellsReflecting( momentumX[step], 0 );
      copyGhostCellsReflecting( momentumY[step], 0 );
      copyGhostCellsReflecting( bathymetry, 20 );
    }
  }

to use the copy function for initialization of the grid like boundary structure: 

.. code-block:: c++

  void WavePropagation2d::copyGhostCellsReflecting( real ** out_grid, real in_value ) {
    idx xMax = cellCountX+1;
    idx yMax = cellCountY+1;

    for( idx x = 1; x < xMax; x++ ) {
      out_grid[x][0] = in_value;
      out_grid[x][yMax] = in_value;
    }

    for( idx y = 1; y < yMax; y++ ) {
      out_grid[0][y] = in_value;
      out_grid[xMax][y] = in_value;
    }

    out_grid[0][0] = in_value;
    out_grid[xMax][0] = in_value;
    out_grid[0][yMax] = in_value;
    out_grid[xMax][yMax] = in_value;
  }

Circular two-dimensional Dambreak
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test the new two-dimensionality of the solver, a circular dam break setup with following inital values: 

 .. math::
  
  \begin{cases} [h, hu, hv]^T = [10, 0, 0]^T &\text{if } \sqrt{x^2+y^2} < 10 \\ 
                [h, hu, hv]^T = [5, 0, 0]^T  \quad &\text{else}
  \end{cases}

is implemented in the computational domain :math:`[-50, 50]^2` :

.. code-block:: c++

  tsunami_lab::setups::DamBreak2d::DamBreak2d( real in_heightInner, real in_heightOuter, real in_radiusDam, real in_xMax, real in_yMax ) {
    heightInner = in_heightInner;
    heightOuter = in_heightOuter;
    radiusDam = in_radiusDam;

    centerDam[0] = in_xMax/2;
    centerDam[1] = in_yMax/2;
  }

  tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getHeight( t_real in_x, t_real in_y ) const {
    real distanceFromCenter = sqrt((centerDam[0]-in_x)*(centerDam[0]-in_x) + (centerDam[1]-in_y)*(centerDam[1]-in_y));
    if( distanceFromCenter < radiusDam ) {
      return heightInner;
    }
    else {
      return heightOuter;
    }
  }

  tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getMomentumX( t_real, t_real ) const {
    return 0;
  }

  tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getMomentumY( t_real, t_real ) const {
    return 0;
  }

  tsunami_lab::t_real tsunami_lab::setups::DamBreak2d::getBathymetry( t_real, t_real ) const {
    return 0;
  }

Visualized solution of 2dDambreak in paraview:

.. video:: _static/dambreak_2d.mp4
  :autoplay:
  :loop:
  :height: 300
  :width: 650

Bathymetry Support in two-dimensional Solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added the bathymetry of the computational domain in the Bathymetry2D setup (added an obstacle):

.. code-block:: c++

  tsunami_lab::setups::Bathymetry2d::Bathymetry2d( real in_heightInner, real in_heightOuter, real in_radiusDam, real in_xMax, real in_yMax, real in_scaling ) {
    heightInner = in_heightInner;
    heightOuter = in_heightOuter;
    radiusDam = in_radiusDam * in_scaling;
    scaling = in_scaling;

    centerDam[0] = (in_xMax/2) * scaling;
    centerDam[1] = (in_yMax/2) * scaling;
  }

  tsunami_lab::t_real tsunami_lab::setups::Bathymetry2d::getHeight( t_real in_x, t_real in_y ) const {
    real distanceFromCenter = sqrt((centerDam[0]-in_x)*(centerDam[0]-in_x) + (centerDam[1]-in_y)*(centerDam[1]-in_y));
    // t_real distanceFromCenter = sqrt(in_x*in_x + in_y*in_y);
    if( distanceFromCenter < radiusDam ) {
      return heightInner;
    }
    else {
      return heightOuter;
    }
  }

  tsunami_lab::t_real tsunami_lab::setups::Bathymetry2d::getMomentumX( t_real, t_real ) const {
    return 0;
  }

  tsunami_lab::t_real tsunami_lab::setups::Bathymetry2d::getMomentumY( t_real, t_real ) const {
    return 0;
  }

  tsunami_lab::t_real tsunami_lab::setups::Bathymetry2d::getBathymetry( t_real in_x, t_real in_y ) const {
    if ( in_x > 7.5 && in_x < 8 && in_y > 7.5 && in_y < 8 ) {
      return -1;
    }
    return -2;
  }

Visualized solution of Bathymetry2D with obstacle (sadly no waves are influenced by bathymetry):

.. video:: _static/bathymetry_bump_2d.mp4
  :autoplay:
  :loop:
  :height: 300
  :width: 650

Individual Member Contributions
--------------------------------

Marek Sommerfeld: two-dimensional solver integration, project report

Moritz Rätz: project report, two-dimensional solver integration


