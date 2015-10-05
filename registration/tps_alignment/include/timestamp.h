#ifndef TIMESTAMP_H
#define TIMESTAMP_H
/*
Copyright 2007 Szymon Rusinkiewicz
Princeton University

timestamp.h
Wrapper around system-specific timestamps.

-----

This file is part of tps_alignment.

tps_alignment is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

tps_alignment is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#ifdef WIN32

# include <limits.h>
# include <windows.h>
# define usleep(x) Sleep((x)/1000)

  struct timestamp { LARGE_INTEGER t; };

  static inline double LI2d(const LARGE_INTEGER &li)
  {
	// Work around random compiler bugs...
	double d = *(unsigned *)(&(li.HighPart));
	d *= 65536.0 * 65536.0;
	d += *(unsigned *)(&(li.LowPart));
	return d;
  }

  static inline float operator - (const timestamp &t1, const timestamp &t2)
  {
	static LARGE_INTEGER PerformanceFrequency;
	static int status = QueryPerformanceFrequency(&PerformanceFrequency);
	if (status == 0) return 1.0f;

	return (LI2d(t1.t) - LI2d(t2.t)) / LI2d(PerformanceFrequency);
  }

  static inline timestamp now()
  {
	timestamp t;
	QueryPerformanceCounter(&t.t);
	return t;
  }

#else

# include <sys/time.h>
# include <unistd.h>

  typedef struct timeval timestamp;

  static inline float operator - (const timestamp &t1, const timestamp &t2)
  {
	return (float)(t1.tv_sec  - t2.tv_sec) +
	       1.0e-6f*(t1.tv_usec - t2.tv_usec);
  }

  static inline timestamp now()
  {
	timestamp t;
	gettimeofday(&t, 0);
	return t;
  }

#endif


#endif
