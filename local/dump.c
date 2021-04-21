pthread_barrier_wait(&bar);
    pthread_mutex_lock(&do_transpose);
    if(!is_tranposed){
        transpose_matrix(v, n);
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
        int i = find_free_bit(rows_bitmap);
        if(i == -1){
            printf("Invalid free row\n");
        }
        rows_bitmap[i] = true;
        pthread_mutex_unlock(&bitmap_mut);
        if(!fft_r2(v+i*n, n, exp)){
            printf("failed\n");
        }
        pthread_mutex_lock(&bitmap_mut);
    }
    pthread_mutex_unlock(&bitmap_mut);


    pthread_barrier_wait(&bar);
    pthread_mutex_lock(&do_transpose2);
    if(!is_tranposed2){
        transpose_matrix(v, n);
        is_tranposed2 = true;
    }
    pthread_mutex_unlock(&do_transpose2);





    #include "fft_lib.h"
#include <pthread.h>
#include <stdint.h>
#ifdef __APPLE__
/* Shim for Mac OS X (use at your own risk ;-) */
# include "apple_pthread_barrier.h"
#endif /* __APPLE__ */
#define NUM_THREADS 3
#define SIZE 4096

#include <string.h>

void transpose_scalar_block(double complex *A, double complex *B, const int lda, const int ldb, const int block_size) {
    for(int i=0; i<block_size; i++) {
        for(int j=0; j<block_size; j++) {
            B[j*ldb + i] = A[i*lda +j];
        }
    }
}

void transpose_block(double complex *A, double complex *B, const int n, const int m, const int lda, const int ldb, const int block_size) {
    for(int i=0; i<n; i+=block_size) {
        for(int j=0; j<m; j+=block_size) {
            transpose_scalar_block(&A[i*lda +j], &B[j*ldb + i], lda, ldb, block_size);
        }
    }
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

	for (i = 0; i < N; i++)
	{
		for (j=i; j < N; j++)
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
        int i = find_free_bit(rows_bitmap);
        if(i == -1){
            printf("Invalid free row\n");
        }
        rows_bitmap[i] = true;
        // printf("\nthread %d doing row %d\n", pthread_self(), i);
        pthread_mutex_unlock(&bitmap_mut);
        if(!fft_r2(v+i*n, n, exp)){
            printf("failed\n");
        }
        pthread_mutex_lock(&bitmap_mut);
    }
    pthread_mutex_unlock(&bitmap_mut);
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

    // double complex *out = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE );
    // transpose_block(vec, out, SIZE, SIZE, SIZE, SIZE, 32);
    // memcpy(vec, out, sizeof(double complex) * SIZE * SIZE);

    transpose_matrix(vec, SIZE);
    for(int i = 0; i < SIZE; ++i){
        rows_bitmap[i] = false;
    }

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
    transpose_matrix(vec, SIZE);
    // memset(out, 0, sizeof(double complex) * SIZE * SIZE);
    // transpose_block(vec, out, SIZE, SIZE, SIZE, SIZE, 32);
    // memcpy(vec, out, sizeof(double complex) * SIZE * SIZE);
    pthread_mutex_destroy(&bitmap_mut);
    pthread_mutex_destroy(&do_transpose);
    pthread_mutex_destroy(&do_transpose2);
    pthread_mutex_destroy(&clear_bitmap);
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
		exptable[i] = cexp(2 * M_PI * i / n * I);

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

    printf("\n");

	printf("Time for FFT was: %.9f\n\n", interval(time_start, time_stop));
    return 0;
}

