#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#define M_PI 3.14159265358979323846


#define SIZE 8192 //Change back to 1<<6
#define SAMPLING_RATE 100
#define FREQUENCY 2


double interval(struct timespec start, struct timespec end)
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


static size_t reverse_bits(size_t val, int width) {
	size_t result = 0;
	int i;
	for (i = 0; i < width; i++, val >>= 1)
		result = (result << 1) | (val & 1U);
	return result;
}


bool Fft_transformRadix2(double complex *vec, size_t n, bool inverse, double complex *exptable) {
	// Length variables
	int levels = 0;  

	size_t temp;
	size_t i;
	size_t j;
	size_t size;
	size_t k;
	size_t m;

	//determine the number of stages in the DIT (log2(n))
	for (temp = n; temp > 1U; temp >>= 1)
		levels++;
	if ((size_t)1U << levels != n)
		return false;  // n is not a power of 2
	
	// // Trigonometric tables
	// if (SIZE_MAX / sizeof(double complex) < n / 2)
	// 	return false;
	
	// Bit-reversed addressing permutation
	// 000 => 000
	// 001 => 100
	// 010 => 010
	// 011 => 110
	// ...
	// for n=8: [a0, a1, a2, a3, a4, a5, a6, a7] => [a0, a4, a2, a6, a1, a5, a2, a7]
	size_t i_index;
	size_t j_index;
	size_t l_index;

	for (m = 0; m < n; m++) {
		for (i = 0; i < n; i++) {
			size_t j = reverse_bits(i, levels);

			if (j > i) {
				i_index = m*n + i;
				j_index = m*n + j;
				double complex temp = vec[i_index];
				vec[i_index] = vec[j_index];
				vec[j_index] = temp;
			}
			// if (j > i) {
			// 	double complex temp = vec[i];
			// 	vec[i] = vec[j];
			// 	vec[j] = temp;
			// }
		}
		
		// Cooley-Tukey decimation-in-time radix-2 FFT
		// loop through each stage
		for (size = 2; size <= n; size *= 2) {														
			size_t halfsize = size / 2;																			
			size_t tablestep = n / size;		

			//for each stage, compute butterly for 2 outputs in groups of 2, 4, 8, ...															
			for (i = 0; i < n; i += size) {	
				// compute butterfly (2 outputs)														

				for (j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
					size_t l = j + halfsize;
					j_index = m*n + j;
					l_index = m*n + l;									
					double complex temp = vec[l_index] * exptable[k];													
					vec[l_index] = vec[j_index] - temp;
					vec[j_index] += temp;
				}

				// for (j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				// 	size_t l = j + halfsize;
				// 	//j_index = m*n + j;
				// 	//l_index = m*n + l;									
				// 	double complex temp = vec[l] * exptable[k];													
				// 	vec[l] = vec[j] - temp;
				// 	vec[j] += temp;
				}
			if (size == n)  // Prevent overflow in 'size *= 2'
				break;
		}
	}
	return true;
}

static void generate_sin_points(double *v, unsigned long size, unsigned long sampling_rate, unsigned long frequency){
    double time = 0.0;
    double inc = (double)1/sampling_rate;
    double W = (double)2 * M_PI * frequency;

    int i;

    for(i = 0; i < size; ++i){
        v[i] = 100*sin(time*W);
        time += inc;
    }
}

static void generate_sin_points_c(double complex *v, unsigned long size, unsigned long sampling_rate, unsigned long frequency){
    double time = 0.0;
    double inc = (double)1/sampling_rate;
    double W = (double)2 * M_PI * frequency;
    int i;

    for(i = 0; i < size; ++i){
        v[i] = (double complex) 100*sin(time*W);
		//printf("v[%d] is %.2lf j%.2lf\n", i, creal(v[i]), cimag(v[i])); //Used for debugging, checked why vector results in function didn't match result outside function
        time += inc;
    }
}
void DFT(double complex *in, double complex *out, unsigned long N){
  int i;
  int j;
  int k;

	// 1D DFT
		// for(i = 0; i < N; ++i){
		// 	double complex sum = 0;
		// 	for(j = 0; j < N; ++j){
		// 		sum += in[j] * cexp((-I*2*M_PI) / N * i * j);
		// 	}
		// 	out[i] = sum;
		// }

	//Version for 2D DFT
	for(k = 0; k < N; ++k){
		for(i = 0; i < N; ++i){
			double complex sum = 0;
			for(j = 0; j < N; ++j){
				sum += in[k*N+j] * cexp((-I*2*M_PI) / N * i * j);
			}
			out[k*N+i] = sum;
		}
	}

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


int main(int argc, char **argv){
    double complex *in = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE);
    double complex *inp = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE );
    double complex *out = (double complex *)malloc(sizeof(double complex) * SIZE * SIZE);

    //generate_sin_points_c(in, SIZE * SIZE, SAMPLING_RATE, FREQUENCY);
    //generate_sin_points_c(inp, SIZE * SIZE, SAMPLING_RATE, FREQUENCY);
    
	printf("SIZE by SIZE is: %d\n", SIZE * SIZE);

	//Generating a matrix of values to do FFT on
    generate_sin_points_c(inp, SIZE * SIZE, SAMPLING_RATE, FREQUENCY);
	
	
    int i, j;
	size_t n = SIZE;

	/*
    for(i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			//printf("i is: %d, j is %d... ", i, j);
        	printf("%.2lf j%.2lf  ", creal(inp[i*n+j]), cimag(inp[i*n+j]));
		}
		printf("\n");
		//printf("%.2lf j%.2lf  ", creal(in[i]), cimag(in[i]));
    }
	*/

    printf("\n");
    printf("\n");
    
	// DFT(in, out, SIZE);

	// for(i = 0; i < SIZE; ++i){
	// 	for (j = 0; j < SIZE; ++j){
    //     	printf("%.2lf j%.2lf  ", creal(out[i*SIZE+j]), cimag(out[i*SIZE+j]));
	// 	}
	// 	printf("\n");
	// 	// printf("%.2lf j%.2lf  ", creal(out[i]), cimag(out[i]));
    // }
    
	// printf("\n");

	bool inverse = 0;

	//set up the table [W^0, W^1, W^2, ..., W^(n/2 - 1)] where W = e^(-2*pi*i/n) for forward fft and W = e^(2*pi*i/n) for inverse fft
	
	double complex *exptable = malloc((SIZE / 2) * sizeof(double complex));
	if (exptable == NULL)
		return false;

	for (i = 0; i < n / 2; i++)
		exptable[i] = cexp((inverse ? 2 : -2) * M_PI * i / n * I);

	// //For 2D, need to do FFT, then transpose, then do FFT again, then transpose back.
    
	printf("Starting now\n");
	struct timespec time_start, time_stop;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
	bool ret = Fft_transformRadix2(inp, SIZE, inverse, exptable);
	transpose_matrix(inp, SIZE);
	ret = Fft_transformRadix2(inp, SIZE, inverse, exptable);
	transpose_matrix(inp, SIZE);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
	printf("Finishing now\n");

    printf("\n");
    
	/*
    if (ret){
      for(i = 0; i < SIZE; ++i){
		  for (j = 0; j < SIZE; ++j){
        	printf("%.2lf j%.2lf   ", creal(inp[i*SIZE+j]), cimag(inp[i*SIZE+j]));
		 }
		 printf("\n");
		//printf("%.2lf j%.2lf   ", creal(inp[i]), cimag(inp[i]));
      }
    }
    else{
      printf("There was an error with your input\n");
    }
	*/

    printf("\n");
	printf("Time for FFT was: %.9f\n\n", interval(time_start, time_stop));
}
