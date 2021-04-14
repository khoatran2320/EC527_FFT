#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#define M_PI 3.14159265358979323846


#define SIZE 1<<6
#define SAMPLING_RATE 100
#define FREQUENCY 2

static size_t reverse_bits(size_t val, int width) {
	size_t result = 0;
	for (int i = 0; i < width; i++, val >>= 1)
		result = (result << 1) | (val & 1U);
	return result;
}


bool Fft_transformRadix2(double complex *vec, size_t n, bool inverse) {
	// Length variables
	int levels = 0;  // Compute levels = floor(log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if ((size_t)1U << levels != n)
		return false;  // n is not a power of 2
	
	// Trigonometric tables
	if (SIZE_MAX / sizeof(double complex) < n / 2)
		return false;
	double complex *exptable = malloc((n / 2) * sizeof(double complex));
	if (exptable == NULL)
		return false;
	for (size_t i = 0; i < n / 2; i++)
		exptable[i] = cexp((inverse ? 2 : -2) * M_PI * i / n * I);
	
	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverse_bits(i, levels);
		if (j > i) {
			double complex temp = vec[i];
			vec[i] = vec[j];
			vec[j] = temp;
		}
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size_t size = 2; size <= n; size *= 2) {														
		size_t halfsize = size / 2;																			
		size_t tablestep = n / size;																	
		for (size_t i = 0; i < n; i += size) {																
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {							
				size_t l = j + halfsize;																	
				double complex temp = vec[l] * exptable[k];													
				vec[l] = vec[j] - temp;
				vec[j] += temp;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
	
	free(exptable);
	return true;
}

static void generate_sin_points(double *v, unsigned long size, unsigned long sampling_rate, unsigned long frequency){
    double time = 0.0;
    double inc = (double)1/sampling_rate;
    double W = (double)2 * M_PI * frequency;

    for(int i = 0; i < size; ++i){
        v[i] = 100*sin(time*W);
        time += inc;
    }
}

static void generate_sin_points_c(double complex *v, unsigned long size, unsigned long sampling_rate, unsigned long frequency){
    double time = 0.0;
    double inc = (double)1/sampling_rate;
    double W = (double)2 * M_PI * frequency;

    for(int i = 0; i < size; ++i){
        v[i] = 100*sin(time*W);
        time += inc;
    }
}
void DFT(double complex *in, double complex *out, unsigned long N){

    for(int i = 0; i < N; ++i){
        double complex sum = 0;
        for(int j = 0; j < N; ++j){
            sum += in[j] * cexp((-I*2*M_PI) / N * i * j);
        }
        out[i] = sum;
    }
}
int main(int argc, char **argv){
    double complex *in = (double complex *)malloc(sizeof(double complex) * SIZE);
    double complex *inp = (double complex *)malloc(sizeof(double complex) * SIZE);
    double complex *out = (double complex *)malloc(sizeof(double complex) * SIZE);
    generate_sin_points_c(in, SIZE, SAMPLING_RATE, FREQUENCY);
    generate_sin_points_c(inp, SIZE, SAMPLING_RATE, FREQUENCY);
    DFT(in, out, SIZE);

    // for(int i = 0; i < SIZE; ++i){
    //     printf("%.2lf j%.2lf   ", creal(in[i]), cimag(in[i]));
    // }
    // printf("\n");
    // printf("\n");
    // for(int i = 0; i < SIZE; ++i){
    //     printf("%.2lf j%.2lf   ", creal(out[i]), cimag(out[i]));
    // }
    // printf("\n");

    bool ret = Fft_transformRadix2(inp, SIZE, 0);
    // printf("\n");
    // for(int i = 0; i < SIZE; ++i){
    //     printf("%.2lf j%.2lf   ", creal(inp[i]), cimag(inp[i]));
    // }
    // printf("\n");
}