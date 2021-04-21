// #define _XOPEN_SOURCE 600
#include "fft_lib.h"
#include <pthread.h>
#include <stdint.h>
#define NUM_THREADS 2
#define SIZE 8192
#include <sys/time.h>
#include "string.h"
#ifdef __APPLE__
/* Shim for Mac OS X (use at your own risk ;-) */
# include "apple_pthread_barrier.h"
#endif /* __APPLE__ */
#include <omp.h>

static inline int rand_int(int lower, int upper)
{
    return rand() % (upper - lower + 1);
}
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
static inline void transpose_scalar_block(double complex *A, double complex *B, const int lda, const int ldb, const int block_size) {
    #pragma omp parallel for
    for(int i=0; i<block_size; i++) {
        for(int j=0; j<block_size; j++) {
            B[j*ldb + i] = A[i*lda +j];
        }
    }
}

static inline void transpose_block(double complex *A, double complex *B, const int n, const int m, const int lda, const int ldb, const int block_size) {
    #pragma omp parallel for
    for(int i=0; i<n; i+=block_size) {
        for(int j=0; j<m; j+=block_size) {
            transpose_scalar_block(&A[i*lda +j], &B[j*ldb + i], lda, ldb, block_size);
        }
    }
}

double wakeup_delay()
{
  double meas = 0; int i, j;
  struct timespec time_start, time_stop;
  double quasi_random = 0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
  j = 100;
  while (meas < 1.0) {
    for (i=1; i<j; i++) {
      /* This iterative calculation uses a chaotic map function, specifically
         the complex quadratic map (as in Julia and Mandelbrot sets), which is
         unpredictable enough to prevent compiler optimisation. */
      quasi_random = quasi_random*quasi_random - 1.923432;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    meas = interval(time_start, time_stop);
    j *= 2; /* Twice as much delay next time, until we've taken 1 second */
  }
  return quasi_random;
}


static int find_free_bit_d(bool *bitmap, int start, int direction){
    if(direction){
        for(int i = start; i < SIZE; ++i){
            if(!bitmap[i]){
                return i;
            }
        }
        for(int i = start - 1; i >= 0; --i){
            if(!bitmap[i]){
                return i;
            }
        }
    }
    else{
        for(int i = start; i >= 0; --i){
            if(!bitmap[i]){
                return i;
            }
        }
        for(int i = start + 1; i < SIZE; ++i){
            if(!bitmap[i]){
                return i;
            }
        }
    }
    return -1;
}

static int find_free_bit(bool *bitmap){
    for(int i = 0; i < SIZE; ++i){
        if(!bitmap[i]){
            return i;
        }
    }

    return -1;
}

static bool is_full(bool *bitmap){
    for(unsigned short i = 0; i < SIZE;++i){
        if(!bitmap[i]){
            return false;
        }
    }
    return true;
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

pthread_mutex_t bitmap_mut;
pthread_mutex_t print_m;
pthread_mutex_t do_transpose;
pthread_mutex_t do_transpose2;
pthread_mutex_t clear_bitmap;


bool is_tranposed = false;
bool is_tranposed2 = false;
bool is_bitmap_cleared = false;
pthread_barrier_t bar;
bool rows_bitmap[SIZE];
struct args{
    double complex *v;
    size_t n;
    double complex *exp;
};


void * start_routine(void *arg){
    struct args *info = (struct args *)arg; 
    double complex *exp = info->exp;
    double complex *v = info->v;
    size_t n = info->n;
    // printf("%ld\n", n);
    // pthread_mutex_lock(&print_m);
    // printf("\n");
    // printf("\n");
    // printf("\n");
    // print_matlab(v, n);
    // pthread_mutex_unlock(&print_m);

    pthread_mutex_lock(&bitmap_mut);
    while(!is_full(rows_bitmap)){
        int i = find_free_bit_d(rows_bitmap, rand_int(0, n-1), rand_int(0,1));
        // int i = find_free_bit(rows_bitmap);
        if(i == -1){
            printf("Invalid free row\n");
            break;
        }
        rows_bitmap[i] = true;
        // printf("\nthread %ld doing row %d\n", pthread_self(), i);
        pthread_mutex_unlock(&bitmap_mut);
        double complex arr[n];
        memcpy(arr, v+i*n, sizeof(double complex)*n);
        if(!fft_r2(arr, n, exp)){
            printf("failed\n");
        }
        memcpy(v+i*n, arr, sizeof(double complex)*n);
        pthread_mutex_lock(&bitmap_mut);
    }
    pthread_mutex_unlock(&bitmap_mut);

    pthread_barrier_wait(&bar);
    pthread_mutex_lock(&do_transpose);
    if(!is_tranposed){
        struct timespec time_start, time_stop;
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


    pthread_mutex_lock(&clear_bitmap);
    if(!is_bitmap_cleared){
        for(int i = 0; i < SIZE; ++i){
            rows_bitmap[i] = false;
        }
        is_bitmap_cleared = true;
    }
    pthread_mutex_unlock(&clear_bitmap);


    pthread_mutex_lock(&bitmap_mut);
    while(!is_full(rows_bitmap)){
        int i = find_free_bit_d(rows_bitmap, rand_int(0, n-1), rand_int(0,1));
        // int i = find_free_bit(rows_bitmap);
        // printf("\n\n\nthread %ld doing row %d\n", pthread_self(), i);
        if(i == -1){
            printf("Invalid free row\n");
            break;
        }
        rows_bitmap[i] = true;
        pthread_mutex_unlock(&bitmap_mut);
        double complex arr[n];
        memcpy(arr, v+i*n, sizeof(double complex)*n);
        if(!fft_r2(arr, n, exp)){
            printf("failed\n");
        }
        memcpy(v+i*n, arr, sizeof(double complex)*n);
        pthread_mutex_lock(&bitmap_mut);
    }
    pthread_mutex_unlock(&bitmap_mut);


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
    pthread_exit(NULL);

}
bool mtfft(double complex *vec, size_t width, double complex *exptable){

    pthread_t tid[NUM_THREADS];
    pthread_mutex_init(&bitmap_mut,  NULL);
    pthread_mutex_init(&do_transpose,  NULL);
    pthread_mutex_init(&do_transpose2,  NULL);
    pthread_mutex_init(&clear_bitmap,  NULL);
    pthread_mutex_init(&print_m,  NULL);
    pthread_barrier_init(&bar, NULL, NUM_THREADS);
    for(int i = 0; i < SIZE; ++i){
        rows_bitmap[i] = false;
    }


    // printf("%ld\n", width);
    struct args arg = (struct args){.exp = exptable, .v = vec, .n = width};
    // printf("%ld\n", arg.n);
    // print_matlab(arg.v, arg.n);

    for(int i = 0; i < NUM_THREADS; ++i){
        if(pthread_create(&tid[i], NULL, start_routine, (void*)&arg)){
            printf("Unable to create thread\n");
            return false;
        }
    }

    for(int i = 0; i < NUM_THREADS; ++i){
        if(pthread_join(tid[i], NULL)){
            printf("Unable to join thread\n");
            return false;
        }
    }
    

    pthread_mutex_destroy(&bitmap_mut);
    pthread_mutex_destroy(&do_transpose);
    pthread_mutex_destroy(&do_transpose2);
    pthread_mutex_destroy(&clear_bitmap);
    pthread_barrier_destroy(&bar);

    return true;
}

int main(int argc, char **argv){

    double temp = wakeup_delay();
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
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);

    clock_t start = clock(), diff;
    bool ret = mtfft(inp, n, exptable);
    diff = clock() - start;
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);


    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds", msec/1000, msec%1000);
    // printf("\n");
    // printf("\n");
    // printf("\n");
    // printf("\n");
    // print_matlab(inp, SIZE);
    // printf("\n");

    // printf("\n");

	// printf("Time for FFT was: %.6lf\n\n", interval(time_start, time_stop));
    return 0;
}

