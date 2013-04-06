#ifdef TITAN

#define WRDLEN WRDLEN

#else
#if defined (AIX) || defined (__hpux)

#define WRDLEN wrdlen

#else

#define WRDLEN wrdlen_

#endif
#endif

#if defined (ULTRIX) || defined (IRIX) || defined (AIX) || defined (LINUX)
#define LWORD  4
#else
#define LWORD  1
#endif


void WRDLEN (int *);

void WRDLEN (int *len)
{
   *len = LWORD;
   return;
}
