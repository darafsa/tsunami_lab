/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * Constants / typedefs used throughout the code.
 **/
#ifndef TSUNAMI_LAB_CONSTANTS_H
#define TSUNAMI_LAB_CONSTANTS_H

#include <cstddef>

namespace tsunami_lab {
  //! integral type for cell-ids, pointer arithmetic, etc.
  typedef std::size_t t_idx;
  typedef t_idx idx;

  //! floating point type
  typedef float t_real;
  typedef t_real real;

  enum Solver { ROE, FWAVE };
  enum Boundary { OUTFLOW, REFLECTING };
}

#endif