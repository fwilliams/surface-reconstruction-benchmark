// ============================================================================

#include <OpenMesh/Core/System/config.hh>
#include <OpenMesh/Tools/Utils/conio.hh>

// ----------------------------------------------------------------- Win32 ----
#ifdef WIN32

#include <conio.h>

namespace OpenMesh {
namespace Utils {

int kbhit()  { return ::kbhit();  }
int getch()  { return ::getch();  }
int getche() { return ::getche(); }

} // Tools
} // AS
// ----------------------------------------------------------------- Other ----
#else

// Based on code published by Floyd Davidson in a newsgroup.

#include <stdio.h>     /* stdout, fflush() */
#if !defined(POSIX_1003_1_2001)
#  include <fcntl.h>
#  include <unistd.h>  
#else
#  include <select.h>  /* select()       */
#endif
#include <termios.h>   /* tcsetattr()    */
#include <sys/ioctl.h> /* ioctl()        */
#include <sys/time.h>  /* struct timeval */

namespace OpenMesh {
namespace Utils {

#ifdef CTIME
#  undef CTIME
#endif
#define CTIME 1
#define CMIN  1


int kbhit(void)
{
  int cnt = 0;
  int error;
  static struct termios Otty, Ntty;

  tcgetattr(0, &Otty);
  Ntty = Otty;

  Ntty.c_iflag      = 0; /* input mode */
  Ntty.c_oflag      = 0; /* output mode */
  Ntty.c_lflag     &= ~ICANON; /* raw mode */
  Ntty.c_cc[VMIN]   = CMIN; /* minimum chars to wait for */
  Ntty.c_cc[VTIME]  = CTIME; /* minimum wait time */

  if (0 == (error = tcsetattr(0, TCSANOW, &Ntty))) 
  {
    struct timeval tv;
    error += ioctl(0, FIONREAD, &cnt);
    error += tcsetattr(0, TCSANOW, &Otty);
    tv.tv_sec = 0;
    tv.tv_usec = 100; /* insert at least a minimal delay */
    select(1, NULL, NULL, NULL, &tv);
  }
  return (error == 0 ? cnt : -1 );
}


int getch(void)
{
  char ch;
  int error;
  static struct termios Otty, Ntty;
  
  fflush(stdout);
  tcgetattr(0, &Otty);
  Ntty = Otty;
  
  Ntty.c_iflag     = 0;        // input mode
  Ntty.c_oflag     = 0;        // output mode
  Ntty.c_lflag    &= ~ICANON;  // line settings  
  Ntty.c_lflag    &= ~ECHO;    // enable echo
  Ntty.c_cc[VMIN]  = CMIN;     // minimum chars to wait for
  Ntty.c_cc[VTIME] = CTIME;    // minimum wait time

  // Conditionals allow compiling with or without flushing pre-existing
  // existing buffered input before blocking.
#if 1
  // use this to flush the input buffer before blocking for new input
#  define FLAG TCSAFLUSH
#else
  // use this to return a char from the current input buffer, or block if
  // no input is waiting.
#  define FLAG TCSANOW
#endif

  if (0 == (error = tcsetattr(0, FLAG, &Ntty))) 
  {
    error = read(0, &ch, 1 );           // get char from stdin
    error += tcsetattr(0, FLAG, &Otty); // restore old settings
  }
  return (error == 1 ? (int) ch : -1 );
}


int getche(void)
{
  char ch;
  int error;
  static struct termios Otty, Ntty;
  
  fflush(stdout);
  tcgetattr(0, &Otty);
  Ntty = Otty;
  
  Ntty.c_iflag     = 0;        // input mode
  Ntty.c_oflag     = 0;        // output mode
  Ntty.c_lflag    &= ~ICANON;  // line settings  
  Ntty.c_lflag    |= ECHO;     // enable echo
  Ntty.c_cc[VMIN]  = CMIN;     // minimum chars to wait for
  Ntty.c_cc[VTIME] = CTIME;    // minimum wait time

  // Conditionals allow compiling with or without flushing pre-existing
  // existing buffered input before blocking.
#if 1
  // use this to flush the input buffer before blocking for new input
#  define FLAG TCSAFLUSH
#else
  // use this to return a char from the current input buffer, or block if
  // no input is waiting.
#  define FLAG TCSANOW
#endif

  if (0 == (error = tcsetattr(0, FLAG, &Ntty))) {
    error = read(0, &ch, 1 );           // get char from stdin
    error += tcsetattr(0, FLAG, &Otty); // restore old settings
  }

  return (error == 1 ? (int) ch : -1 );
}

} // namespace Tools
} // namespace AS
// ----------------------------------------------------------------------------
#endif // System dependent parts
// ============================================================================

//#define Test
#if defined(Test)

#include <ctype.h>

int main (void) 
{ 
   char  msg[] = "press key to continue...";
   char *ptr   = msg, tmp; 

  while ( !OpenMesh::Utils::kbhit() )
  {
    tmp = *ptr;
    *ptr = islower(tmp) ? toupper(tmp) : tolower(tmp);
    printf("\r%s", msg); fflush(stdout);
    *ptr = (char)tmp;
    if (!*(++ptr)) 
      ptr = msg; 
    usleep(20000);
  }

  printf("\r%s.", msg); fflush(stdout);
  OpenMesh::Utils::getch();
  printf("\r%s..", msg); fflush(stdout);
  OpenMesh::Utils::getche();
  return 0;
}

#endif // Test

// ============================================================================
