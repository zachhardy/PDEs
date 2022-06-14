#include "timer.h"

#include <cmath>
#include <ctime>


void Timer::start()
{ start_time = std::chrono::steady_clock::now(); }


void Timer::stop()
{ end_time = std::chrono::steady_clock::now(); }


double
Timer::get_time()
{
  using namespace std::chrono;
  return duration<double, std::milli>(end_time - start_time).count();
}


std::string Timer::get_time_string()
{
  double dt = get_time()/1000.0;
  int hours = std::floor(dt/3600.0);
  int minutes = std::floor((dt - 3600.0*hours)/60.0);
  int seconds = (int)dt - 3600*hours - 60*minutes;

  char buff[100];
  sprintf(buff, "%02d:%02d:%02d", hours, minutes, seconds);
  return std::string(buff);

}