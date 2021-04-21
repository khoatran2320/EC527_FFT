#include "fft_lib.h"
#define SIZE 512
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

int main(int argc, char **argv){
    double complex *in = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE);
    double complex *inp = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE );
    double complex *out = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE);

	struct timespec time_start, time_stop;

    
	printf("SIZE by SIZE is: %d\n", SIZE * SIZE);

	//Generating a matrix of values to do FFT on
    generate_sin_points_c(inp, SIZE * SIZE, SAMPLING_RATE, FREQUENCY);
	
    int i, j;
	size_t n = SIZE;

    // for(i = 0; i < n; i++){
	// 	for (j = 0; j < n; j++){
    //     	printf("%.2lf j%.2lf  ", creal(inp[i*n+j]), cimag(inp[i*n+j]));
	// 	}
	// 	printf("\n");
    // }

    // printf("\n");
    // printf("\n");
    
    // print_matlab(inp, SIZE);

	//set up the table [W^0, W^1, W^2, ..., W^(n/2 - 1)] where W = e^(-2*pi*i/n) for forward fft and W = e^(2*pi*i/n) for inverse fft
	double complex *exptable = malloc((SIZE / 2) * sizeof(double complex));
	if (exptable == NULL)
		return false;

	for (i = 0; i < n / 2; i++)
		exptable[i] = cexp(2 * M_PI * i / n * I);

	// //For 2D, need to do FFT, then transpose, then do FFT again, then transpose back.
    
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
    bool ret = fft2d(inp, n, exptable);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    // print_matlab(inp, SIZE);
    // printf("\n");
    
    // if(ret){
    //     for(i = 0; i < SIZE; ++i){
    //         for (j = 0; j < SIZE; ++j){
    //             printf("%.2lf j%.2lf   ", creal(inp[i*SIZE+j]), cimag(inp[i*SIZE+j]));
    //         }
    //         printf("\n");
    //     }
    // }
    // else{
    //   printf("There was an error with your input\n");
    // }

    printf("\n");
	printf("Time for FFT was: %.9f\n\n", interval(time_start, time_stop));
}