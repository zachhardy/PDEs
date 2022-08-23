#include "timer.h"

#include <cmath>
#include <ctime>


using namespace PDEs;


void
Timer::start()
{
  start_time = std::chrono::steady_clock::now();
}


void
Timer::stop()
{
  end_time = std::chrono::steady_clock::now();
}


double
Timer::get_time()
{
  using namespace std::chrono;
  return duration<double, std::milli>(end_time - start_time).count();
}