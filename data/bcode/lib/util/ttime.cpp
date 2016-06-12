#include "util.h"
  

Timer::Timer()
{
  times(&colin_time);
  tmu = colin_time.tms_utime;
  tms = colin_time.tms_stime;
}

void Timer::start ()
{
  times(&colin_time);
  tmu = colin_time.tms_utime;
  tms = colin_time.tms_stime;
}

void Timer::stop ()
{
  float local_time,utime,stime;
  times (&colin_time);
  utime=(float)(colin_time.tms_utime-tmu);
  stime=(float)(colin_time.tms_stime-tms);
  
  local_time = (float)(utime + stime);
  cout << "total time: " << (local_time /60.0)<< " user "<< (utime/60.0)<< 
    " system "  << (stime/60.0) <<  endl;
}

float Timer::getTime ()
{
  float local_time;
  times (&colin_time);
  local_time = (float)(colin_time.tms_utime + colin_time.tms_stime-tmu-tms);
  return (local_time /60.0);
}
