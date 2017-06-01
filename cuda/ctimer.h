#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

/**
 * Lean performance timer structure.
 */
typedef struct CPerfTimer {
  struct timeval tv;
  unsigned long us;
  const char *tag;
} CPerfTimer;

/**
 * Initialise and return a performance timer with a given tag.
 *
 * @param s tag string to be attached to the timer
 * @return new CPerfTimer
 */
static inline CPerfTimer startTimer(const char *s) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long us = 1000000 * tv.tv_sec + tv.tv_usec;
  CPerfTimer pf = { .tv = tv, .us = us, .tag = s};
  return pf;
}

/**
 * Print time since creation of pt, together with pt's tag. Only prints
 * results if compiled with -DTIMER flag.
 *
 * @param pt CPerfTimer (usually initialised by startTimer())
 */
static inline void stopTimer(CPerfTimer pt) {
#ifdef TIMER
  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long us = 1000000 * tv.tv_sec + tv.tv_usec;
  if (strlen(pt.tag) < 50) {
    char buf[50] = "                                                 ";
    memcpy(buf, pt.tag, strlen(pt.tag));
    printf("%s: %.3fms\n", buf, (us - pt.us)/((double) 1000.0));
  } else {
    printf("%s: %.3fms\n", pt.tag, (us - pt.us)/((double) 1000.0));
  }
#endif
  return;
}

