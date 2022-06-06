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
  Timer() noexcept = default;

  void start();
  void stop();

  double get_time();
  std::string get_time_string();
};


#endif //TIMER_H
