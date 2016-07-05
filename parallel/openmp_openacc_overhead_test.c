#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>

// Define constants for work loop (type "a*x+b")
#define ACONST 2.3f
#define BCONST 7.13f

/*

  Single core version

*/
long compute_singlecore(float * work, long nelements) {

  // Initialise timer
  struct timespec start_time, end_time;
  clock_gettime(CLOCK_REALTIME, &start_time);

  // Do work
  for (long i = 0; i < nelements; i++) {
    work[i] = ACONST*work[i] + BCONST;
  }

  // Read timer again and compute time spent on loop
  clock_gettime(CLOCK_REALTIME, &end_time);
  long time_elapsed_nanos = end_time.tv_nsec - start_time.tv_nsec;;

  // Set huge number if counter returns erroneous results
  if (time_elapsed_nanos < 0) {
    time_elapsed_nanos = LONG_MAX;
  }

  return time_elapsed_nanos;

}

/*

  Multicore version

*/
long compute_multicore(float * work, long nelements) {

  struct timespec start_time, end_time;
  clock_gettime(CLOCK_REALTIME, &start_time);

  // Do work in parallel on multicore CPU
#pragma omp parallel for num_threads(4) schedule(static) \
  shared(work, nelements) default(none)
  for (long i = 0; i < nelements; i++) {
    work[i] = ACONST*work[i] + BCONST;
  }

  clock_gettime(CLOCK_REALTIME, &end_time);
  long time_elapsed_nanos = end_time.tv_nsec - start_time.tv_nsec;;
  if (time_elapsed_nanos < 0) {
    time_elapsed_nanos = LONG_MAX;
  }

  return time_elapsed_nanos;

}

/*

  GPU version

*/
long compute_gpu(float * work, long nelements) {

  // Copy data to GPU memory (not measured by timer)
#pragma acc enter data copyin(work[0:nelements])

  struct timespec start_time, end_time;
  clock_gettime(CLOCK_REALTIME, &start_time);

  // Do work on GPU
#pragma acc kernels loop independent present(work)
  for (long i = 0; i < nelements; i++) {
    work[i] = ACONST*work[i] + BCONST;
  }

  clock_gettime(CLOCK_REALTIME, &end_time);
  long time_elapsed_nanos = end_time.tv_nsec - start_time.tv_nsec;;
  if (time_elapsed_nanos < 0) {
    time_elapsed_nanos = LONG_MAX;
  }

  // Copy data back to host memory (not measured by timer)
#pragma acc exit data copyout(work[0:nelements])

  return time_elapsed_nanos;

}

/*

  Verify results

*/
long checkresult(float * data, float * result, long nelements) {
  long nerr = 0;
  for (long i = 0; i < nelements; i++) {
    if (result[i] != (ACONST*data[i] + BCONST)) {
      nerr++;
    }
  }
  return nerr;
}

int main() {

  /*
    Parameters
  */

  // Number of array size steps
  int size_steps = 70;

  // Number of times a performance measurement is repeated
  // to get statistics
  int nrepeat = 30;

  // Increase work size in loop
  for (int size_step = 1; size_step <= size_steps; size_step++) {

    long ntickssingle = LONG_MAX;
    long nticksmulti = LONG_MAX;
    long nticksgpu = LONG_MAX;

    // Compute decade, start with 10
    long decade = pow(10, (size_step+8)/9);

    // Compute array size
    long nelements = (((size_step+8)%9)+1)*decade;

    // Get memory
    float * randarray = (float*)malloc(nelements*sizeof(float));
    float * workarray = (float*)malloc(nelements*sizeof(float));

    // Fill array with random numbers
    for (long i = 0; i < nelements; i++) {
      randarray[i] = (float)rand();
    }

    /*

      Single core

    */
    for (int iter = 0; iter < nrepeat; iter ++) {

      // Reset work array
      for (long i = 0; i < nelements; i++) {
        workarray[i] = randarray[i];
      }

      // Run computation and check results
      long ntickssingle_iter = compute_singlecore(workarray, nelements);
      if( checkresult(randarray, workarray, nelements) ) {
        fprintf(stderr, "ERROR in single core computation!\n");
      }

      // Get minimum
      if (ntickssingle_iter < ntickssingle) {
        ntickssingle = ntickssingle_iter;
      }

    }
    /*

      Multicore

    */
    for (int iter = 0; iter < nrepeat; iter ++) {

      for (long i = 0; i < nelements; i++) {
        workarray[i] = randarray[i];
      }

      long nticksmulti_iter = compute_multicore(workarray, nelements);
      if( checkresult(randarray, workarray, nelements) ) {
        fprintf(stderr, "ERROR in multicore computation!\n");
      }

      if (nticksmulti_iter < nticksmulti) {
        nticksmulti = nticksmulti_iter;
      }

    }

    /*

      GPU

    */
    for (int iter = 0; iter < nrepeat; iter ++) {

      for (long i = 0; i < nelements; i++) {
        workarray[i] = randarray[i];
      }

      long nticksgpu_iter = compute_gpu(workarray, nelements);
      if( checkresult(randarray, workarray, nelements) ) {
        fprintf(stderr, "ERROR in GPU computation!\n");
      }

      if (nticksgpu_iter < nticksgpu) {
        nticksgpu = nticksgpu_iter;
      }

    }

    // Print results and tidy up
    printf("%i %i %f %f\n", size_step, nelements, ((float)ntickssingle)/((float)nticksmulti),
           ((float)ntickssingle)/((float)nticksgpu));

    free(workarray);
    free(randarray);
  }

}
