#ifndef MACROS_H
#define MACROS_H

#include <cassert>
#include <stdexcept>
#include <sstream>
#include <iostream>

/**
 * Macros to open and close the code base namespace.
 */
#define PDES_NAMESPACE_OPEN namespace PDEs {
#define PDES_NAMESPACE_CLOSE }


/**
 * General assert statement.
 */
#define Assert(cond, exc)                                   \
  {                                                         \
    if (!(cond))                                            \
    {                                                       \
      std::cerr                                             \
        << std::endl                                        \
        << "----------------------------------------"       \
        << std::endl                                        \
        << "Assertion falied at " << __FILE__ << ":"        \
        << __LINE__ << " inside " << __PRETTY_FUNCTION__    \
        << std::endl << "Condition: " << #cond              \
        << std::endl << "Message: " << #exc                 \
        << std::endl                                        \
        << "----------------------------------------"       \
        << std::endl;                                       \
      exit(EXIT_FAILURE);                                   \
    }                                                       \
  }
#endif //MACROS_H
