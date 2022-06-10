#ifndef UTILITIES_H
#define UTILITIES_H

#include <cstdint>
#include <cstddef>


namespace numbers
{
  /**
   * Representation of the largest number that can be stored by an unsigned
   * integer. This is used when necessary to designate an invalid value.
   */
  const unsigned int invalid_unsigned_int = static_cast<unsigned int>(-1);


  /**
   * Representation of the largest number that can be stored by a size_t
   * integer. This is used when necessary to designate an invalid value.
   */
  const size_t invalid_size_t = static_cast<size_t>(-1);
}

#endif //UTILITIES_H
