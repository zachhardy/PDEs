#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <chrono>


/**
 * An object for profiling bits of code.
 */
class Timer
{
private:
  std::chrono::steady_clock::time_point start_time;
  std::chrono::steady_clock::time_point end_time;

public:
  /**
   * Default constructor.
   */
  Timer() noexcept = default;

  /**
   * Start the timer.
   */
  void
  start();

  /**
   * Stop the timer.
   */
  void
  stop();

  /**
   * Get the duration of time between the start and stop calls.
   */
  double
  get_time();

  /**
   * Print the elapsed time between the start and stop calls.
   *
   * \see Timer::get_time
   */
  std::string
  get_time_string();
};


#endif //TIMER_H
