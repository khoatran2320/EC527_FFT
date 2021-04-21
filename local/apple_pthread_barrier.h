/* Originally inspired by this implementation:
     Albert Armea, "Using pthread_barrier on Mac OS X"
     blog.albertarmea.com/post/47089939939/using-pthreadbarrier-on-mac-os-x

   This version has been greatly edited to fix bugs (more error
   checking, always setting errno, etc.) and to add documentation
   (comments). */

/* To view all predefined macros (things defined by a #define) you can
   type the command:

     gcc -dM $OPTIONS -E - < /dev/null | sort

Apple's compilers (installed as part of Xcode) all define __APPLE__
and __APPLE_CC__ and if targeting Mac OS X they also define __MACH__.
I'm using two of these to "protect" this file and prevent it from being
invoked in other systems. (For example, a non-Apple Mach-based system
probably has the complete pthreads implementation)

*/

#ifdef __MACH__
#ifdef __APPLE_CC__

#ifndef _APPLE_PTHREAD_BARRIER_SHIM_
#define _APPLE_PTHREAD_BARRIER_SHIM_

#include <pthread.h>
#include <errno.h>

/* The value of PTHREAD_BARRIER_SERIAL_THREAD that we get on Centos 7 systems
   is -1 */
#define PTHREAD_BARRIER_SERIAL_THREAD -1

typedef int pthread_barrierattr_t;
typedef struct pthread_barrier_t {
  pthread_mutex_t mutex;
  pthread_cond_t condition;
  int got_count;
  int threshold;
} pthread_barrier_t;

/* 

  pthread_barrier_init    - Initialise a barrier object

int pthread_barrier_init(pthread_barrier_t *restrict barrier,
       const pthread_barrierattr_t *restrict attr, unsigned count);

Allocates the needed resources for a "barrier" object (for use by
pthread_barrier_wait(), see below) and sets the threshold (number of needed
waits before all waits get unblocked) to a number given by "count".

For full specification see:
pubs.opengroup.org/onlinepubs/009604599/functions/pthread_barrier_init.html

*/
int pthread_barrier_init(
  pthread_barrier_t *restrict barrier,
  const pthread_barrierattr_t *restrict attr,
  unsigned count)
{
  if (count == 0) {
    /* The threshold must be nonzero */
    errno = EINVAL;
    return -1;
  }
  /* Create the mutex we'll use to manage the counts */
  if (pthread_mutex_init(&barrier->mutex, 0) < 0) {
    /* Couldn't create it, let's assume it's because "out of resources" */
    errno = EAGAIN;
    return -1;
  }
  if (pthread_cond_init(&barrier->condition, 0) < 0) {
    /* Couldn't create it, let's assume it's because "out of resources".
       must destroy the mutex that was just successfully created */
    pthread_mutex_destroy(&barrier->mutex);
    errno = EAGAIN;
    return -1;
  }
  /* Remember how many simultaneous pthread_barrier_wait calls we need
     before letting them all proceed */
  barrier->threshold = count;
  /* We start with zero threads waiting. */
  barrier->got_count = 0;

  return 0;
}

/*
   pthread_barrier_wait   - Synchronize at a barrier

One or more threads (called "participating threads") can call this
function when they want to wait for synchronisation to occur. 

The calling thread shall block until the required number of threads
(equal to the "count" parameter that was given in
pthread_barrier_create) have called pthread_barrier_wait() specifying
the barrier. When that happens, all of the blocked threads get
unblocked, i.e. their calls to pthread_barrier_wait() return.

The return value indicates if the calling thread was "last": the
thread that caused the threshold to be reached gets the return value
PTHREAD_BARRIER_SERIAL_THREAD, while all other threads get 0.

For full specification see:
pubs.opengroup.org/onlinepubs/009604599/functions/pthread_barrier_wait.html

 */
int pthread_barrier_wait(pthread_barrier_t *barrier)
{
  /* Get the mutex (which protects the count) */
  pthread_mutex_lock(&barrier->mutex);
  (barrier->got_count) += 1;
  if (barrier->got_count >= barrier->threshold) {
    barrier->got_count = 0;
    /* Condition has been met, tell everyone */
    pthread_cond_broadcast(&barrier->condition);
    pthread_mutex_unlock(&barrier->mutex);
    return PTHREAD_BARRIER_SERIAL_THREAD;
  } else {
    /* We don't have enough yet, we need to wait until enough other
       threads come along, and eventually one of them will discover the
       condition has been met, and then they'll tell us. */
    pthread_cond_wait(&barrier->condition, &(barrier->mutex))
    /* . . . a timeless interval passes here . . . */
                                                        ; /* It happened! */
    pthread_mutex_unlock(&barrier->mutex);
    return 0;
  }
}

/* 

  pthread_barrier_destroy    - Destroy a barrier object

int pthread_barrier_destroy(pthread_barrier_t *restrict barrier);

If the barrier object is "idle" (no threads waiting), this frees the
resources used by the barrier object. Otherwise, it tries to report an
appropriate error code.

For full specification see:
pubs.opengroup.org/onlinepubs/009604599/functions/pthread_barrier_destroy.html

*/
int pthread_barrier_destroy(pthread_barrier_t *barrier)
{
  int cur_count;
  /* See if there is a barrier wait partially in progress... */
  pthread_mutex_lock(&barrier->mutex);
  if (barrier->got_count) {
    /* Yes, one or more tasks are waiting, the barrier is "busy" */
    pthread_mutex_unlock(&barrier->mutex);
    errno = EBUSY;
    return -1;
  } else {
    /* It's free, so we can tear it down */
    pthread_cond_destroy(&barrier->condition);
    /* We have to unlock the mutex before destroying it. There's a race
       condition here if some thread tries to do a wait before we're
       done destroying it. However this falls under the provision in the
       official specification, "The effect of subsequent use of the
       barrier is undefined until the barrier is reinitialized by
       another call to pthread_barrier_init()." */
    pthread_mutex_unlock(&barrier->mutex);
    pthread_mutex_destroy(&barrier->mutex);
  }
  return 0;
}

#endif // _APPLE_PTHREAD_BARRIER_SHIM_

#endif // __APPLE_CC__
#endif // __MACH__