#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <cassert>
#include <stdexcept>
#include <sstream>
#include <iostream>

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
#endif //EXCEPTIONS_H
