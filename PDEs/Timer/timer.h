#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <chrono>


namespace PDEs
{
  /**
   * An object for profiling bits of code.
   */
  class Timer
  {
  public:
    /**
     * Default constructor.
     */
    Timer() noexcept = default;

    /**
     * Start the timer. This routine sets the \p start_time.
     */
    void
    start();

    /**
     * Stop the timer. This routine sets the \p stop_time
     */
    void
    stop();

    /**
     * Get the amount of time that elapsed between the \ref start and \ref
     * stop calls.
     */
    double
    get_time();

  private:
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point end_time;
  };
}

#endif //TIMER_H
