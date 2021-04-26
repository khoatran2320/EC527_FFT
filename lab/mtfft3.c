// #define _XOPEN_SOURCE 600
#include "fft_lib.h"
#include <pthread.h>
#include <stdint.h>
#define NUM_THREADS 8
#define SIZE 16384
#include <sys/time.h>
#include "string.h"
// #ifdef __APPLE__
// /* Shim for Mac OS X (use at your own risk ;-) */
// # include "apple_pthread_barrier.h"
// #endif /* __APPLE__ */

static double interval(struct timespec start, struct timespec end)
{
  struct timespec temp;
  temp.tv_sec = end.tv_sec - start.tv_sec;
  temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  if (temp.tv_nsec < 0) {
    temp.tv_sec = temp.tv_sec - 1;
    temp.tv_nsec = temp.tv_nsec + 1000000000;
  }
  return (((double)temp.tv_sec) + ((double)temp.tv_nsec)*1.0e-9);
}

void transpose_matrix(double complex *in, unsigned long N){
	int i;
	int j;
	double complex temp;
    // #pragma omp parallel for private(temp, N) shared(in, i, j)
	for (i = 0; i < N-1; i++)
	{
		for (j=i+1; j < N; j++)
		{
			temp = in[i*N+j];
			in[i*N+j] = in[j*N+i];
			in[j*N+i] = temp;
		}
	}
}

bool fft2d(double complex *vec, size_t width, double complex *exptable){
    for(int i = 0; i < width; ++i){
        if(!fft_r2(vec+i*width, width, exptable)){
            return false;
        }
    }
    transpose_matrix(vec, width);
    for(int i = 0; i < width; ++i){
        if(!fft_r2(vec+i*width, width, exptable)){
            return false;
        }
    }
    transpose_matrix(vec, width);
    return true;
}

pthread_mutex_t do_transpose;


bool is_tranposed = false;
pthread_barrier_t bar;
struct args{
    int row_start;
    int row_end;
    double complex *v;
    size_t n;
    double complex *exp;
};


void * start_routine(void *arg){
    struct args *info = (struct args *)arg; 
    double complex *exp = info->exp;
    double complex *v = info->v;
    size_t n = info->n;
    int row_start = info->row_start;
    int row_end = info->row_end;
    
    //1D fft on rows
    for(int i = 0; i < row_end-row_start+1; ++i){
        if(!fft_r2(v+row_start*n+i*n, n, exp)){
            printf("failed\n");
        }
    }
    //wait for all threads to finish
    pthread_barrier_wait(&bar);

    //one thread tranposes the matrix
    pthread_mutex_lock(&do_transpose);
    if(!is_tranposed){
        transpose_matrix(v, n);
        is_tranposed = true;
    }
    pthread_mutex_unlock(&do_transpose);

    //1D fft on rows
    for(int i = 0; i < row_end-row_start+1; ++i){
        if(!fft_r2(v+row_start*n+i*n, n, exp)){
            printf("failed\n");
        }
    }
    pthread_exit(NULL);
}

int main(int argc, char **argv){
    double complex *inp = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE );
    struct args arg[NUM_THREADS];
	struct timespec time_start, time_stop;
    pthread_t tid[NUM_THREADS];

    pthread_mutex_init(&do_transpose,  NULL);
    pthread_barrier_init(&bar, NULL, NUM_THREADS);
	printf("SIZE by SIZE is: %d\n", SIZE * SIZE);

	//Generating a matrix of values to do FFT on
    generate_sin_points_c(inp, SIZE * SIZE, SAMPLING_RATE, FREQUENCY);
	
    int i, j;
	size_t n = SIZE;

	//set up the table [W^0, W^1, W^2, ..., W^(n/2 - 1)] where W = e^(-2*pi*i/n) for forward fft and W = e^(2*pi*i/n) for inverse fft
	double complex *exptable = malloc((SIZE / 2) * sizeof(double complex));
	if (exptable == NULL)
		return false;

	for (i = 0; i < n / 2; i++)
		exptable[i] = cexp(-2 * M_PI * i / n * I);

    //launch threads
	for(int i = 0; i < NUM_THREADS; ++i){
        arg[i] = (struct args){.exp = exptable, .v = inp, .n = SIZE, .row_start = i*SIZE/NUM_THREADS, .row_end = (i+1)*SIZE/NUM_THREADS-1};
        if(pthread_create(&tid[i], NULL, start_routine, (void*)&arg[i])){
            printf("Unable to create thread\n");
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &time_start);

    //allow threads to finish
    for(int i = 0; i < NUM_THREADS; ++i){
        if(pthread_join(tid[i], NULL)){
            printf("Unable to join thread\n");
            return false;
        }
    }

    //do a final tranpose
    transpose_matrix(inp, SIZE);
    clock_gettime(CLOCK_MONOTONIC_RAW, &time_stop);

	printf("Time for FFT was: %.9f\n\n", interval(time_start, time_stop));
    pthread_mutex_destroy(&do_transpose);
    pthread_barrier_destroy(&bar);
    free(exptable);
    return 0;
}

