#define _XOPEN_SOURCE 600
#include "fft_lib.h"
#include <pthread.h>
#include <stdint.h>
#define NUM_THREADS 8
#define SIZE 8192
#include <sys/time.h>
#include "string.h"
#ifdef __APPLE__
/* Shim for Mac OS X (use at your own risk ;-) */
# include "apple_pthread_barrier.h"
#endif /* __APPLE__ */
#include <omp.h>

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
pthread_mutex_t do_transpose2;


bool is_tranposed = false;
bool is_tranposed2 = false;
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
    
    double complex *vec1 = (double complex *)malloc((row_end-row_start)*n*sizeof(double complex));
    memcpy(vec1, v+row_start*n, (row_end-row_start)*n*sizeof(double complex));
    for(int i = 0; i < row_end-row_start; ++i){
        if(!fft_r2(vec1+i*n, n, exp)){
            printf("failed\n");
        }
    }
    memcpy(v+row_start*n, vec1, (row_end-row_start)*n*sizeof(double complex));
    pthread_barrier_wait(&bar);
    pthread_mutex_lock(&do_transpose);
    if(!is_tranposed){
        // struct timespec time_start, time_stop;
        // clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
        // double complex *temp = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE );
        transpose_matrix(v, n);
        // transpose_block(v, temp, n, n, n, n, 4);
        // memcpy(v, temp, sizeof(double complex) * SIZE * SIZE );
        // clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
        // printf("Time for tranpose was: %.6f\n\n", interval(time_start, time_stop));
        is_tranposed = true;
    }
    pthread_mutex_unlock(&do_transpose);

    double complex *vec2 = (double complex *)malloc((row_end-row_start)*n*sizeof(double complex));
    memcpy(vec2, v+row_start*n, (row_end-row_start)*n*sizeof(double complex));
    for(int i = 0; i < row_end-row_start; ++i){
        if(!fft_r2(vec2+i*n, n, exp)){
            printf("failed\n");
        }
    }
    memcpy(v+row_start*n, vec2, (row_end-row_start)*n*sizeof(double complex));

    pthread_barrier_wait(&bar);
    pthread_mutex_lock(&do_transpose2);
    if(!is_tranposed2){
        // double complex *temp = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE );
        transpose_matrix(v, n);
        // transpose_block(v, temp, n, n, n, n, 4);
        // memcpy(v, temp, sizeof(double complex) * SIZE * SIZE );
        is_tranposed2 = true;
    }
    pthread_mutex_unlock(&do_transpose2);

    // free(vec1);
    // free(vec2);
    pthread_exit(NULL);

}
bool mtfft(double complex *vec, size_t width, double complex *exptable){

    pthread_t tid[NUM_THREADS];
    pthread_mutex_init(&do_transpose,  NULL);
    pthread_mutex_init(&do_transpose2,  NULL);
    pthread_barrier_init(&bar, NULL, NUM_THREADS);

    double complex *exptable1 = malloc((SIZE / 2) * sizeof(double complex));
    double complex *exptable2 = malloc((SIZE / 2) * sizeof(double complex));
    double complex *exptable3 = malloc((SIZE / 2) * sizeof(double complex));
    memcpy(exptable1, exptable, (SIZE / 2) * sizeof(double complex));
    memcpy(exptable2, exptable, (SIZE / 2) * sizeof(double complex));
    memcpy(exptable3, exptable, (SIZE / 2) * sizeof(double complex));
    // printf("%ld\n", width);
    struct args arg1 = (struct args){.exp = exptable1, .v = vec, .n = width, .row_start = 0, .row_end = SIZE/3};
    struct args arg2 = (struct args){.exp = exptable2, .v = vec, .n = width, .row_start = SIZE/3, .row_end = 2*SIZE/NUM_THREADS};
    struct args arg3 = (struct args){.exp = exptable3, .v = vec, .n = width, .row_start = 2*SIZE/NUM_THREADS, .row_end = SIZE};
    // printf("%ld\n", arg.n);
    // print_matlab(arg.v, arg.n);

    int i = 0;

    if(pthread_create(&tid[i++], NULL, start_routine, (void*)&arg3)){
        printf("Unable to create thread\n");
        return false;
    }
    if(pthread_create(&tid[i++], NULL, start_routine, (void*)&arg1)){
        printf("Unable to create thread\n");
        return false;
    }
    if(pthread_create(&tid[i++], NULL, start_routine, (void*)&arg2)){
        printf("Unable to create thread\n");
        return false;
    }


    for(int i = 0; i < NUM_THREADS; ++i){
        if(pthread_join(tid[i], NULL)){
            printf("Unable to join thread\n");
            return false;
        }
    }

    free(exptable1);
    free(exptable2);
    free(exptable3);
    pthread_mutex_destroy(&do_transpose);
    pthread_mutex_destroy(&do_transpose2);
    pthread_barrier_destroy(&bar);

    return true;
}

int main(int argc, char **argv){
    double complex *inp = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE );

	struct timespec time_start, time_stop;

    
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

	// //For 2D, need to do FFT, then transpose, then do FFT again, then transpose back.
    
	// print_matlab(inp, SIZE);
    // printf("\n");
    // printf("\n");
    // printf("\n");
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
    bool ret = mtfft(inp, n, exptable);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);

    // printf("\n");
    // printf("\n");
    // printf("\n");
    // printf("\n");
    // print_matlab(inp, SIZE);
    // printf("\n");

    // printf("\n");

	printf("Time for FFT was: %.9f\n\n", interval(time_start, time_stop));
    return 0;
}

