#ifndef FFT_LIB_H
#define FFT_LIB_H
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>


// #define SIZE (1* 1<<10)
#define SAMPLING_RATE 100
#define FREQUENCY 2
#define M_PI 3.14159265358979323846


static inline size_t reverse_bits(size_t val, int width) {
	size_t result = 0;
	for (int i = 0; i < width; i++, val >>= 1)
		result = (result << 1) | (val & 1U);
	return result;
}
bool fft_r2(double complex *vec, size_t n, double complex *exptable) {
	// Length variables
	int levels = 0;  

	//determine the number of stages in the DIT (log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if ((size_t)1U << levels != n){
		return false;  // n is not a power of 2
    }
	
	
	
	// Bit-reversed addressing permutation
	// 000 => 000
	// 001 => 100
	// 010 => 010
	// 011 => 110
	// ...
	// for n=8: [a0, a1, a2, a3, a4, a5, a6, a7] => [a0, a4, a2, a6, a1, a5, a2, a7]
	for (size_t i = 0; i < n; i++) {
		size_t j = reverse_bits(i, levels);
		if (j > i) {
			double complex temp = vec[i];
			vec[i] = vec[j];
			vec[j] = temp;
		}
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	// loop through each stage
	for (size_t size = 2; size <= n; size *= 2) {														
		size_t halfsize = size / 2;																			
		size_t tablestep = n / size;		

		//for each stage, compute butterly for 2 outputs in groups of 2, 4, 8, ...															
		for (size_t i = 0; i < n; i += size) {	
				
			// compute butterfly (2 outputs)														
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
	
	return true;
}

void generate_sin_points(double *v, unsigned long size, unsigned long sampling_rate, unsigned long frequency){
    double time = 0.0;
    double inc = (double)1/sampling_rate;
    double W = (double)2 * M_PI * frequency;

    for(int i = 0; i < size; ++i){
        v[i] = 100*sin(time*W);
        time += inc;
    }
}

void generate_sin_points_c(double complex *v, unsigned long size, unsigned long sampling_rate, unsigned long frequency){
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

void print_matlab(double complex *vec, size_t width){
    printf("[");
    for(int i = 0; i < width;++i){
        for(int j = 0; j < width; ++j){
            printf("%lf ", creal(vec[i*width+j]));
        }
        printf(";\n");
    }
    printf("];");
}

#endif