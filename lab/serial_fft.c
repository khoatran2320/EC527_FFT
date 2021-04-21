#include "fft_lib.h"




int main(int argc, char **argv){
    double complex *in = (double complex *)malloc(sizeof(double complex) * SIZE);
    double complex *inp = (double complex *)malloc(sizeof(double complex) * SIZE);
    double complex *out = (double complex *)malloc(sizeof(double complex) * SIZE);
    generate_sin_points_c(in, SIZE, SAMPLING_RATE, FREQUENCY);
    generate_sin_points_c(inp, SIZE, SAMPLING_RATE, FREQUENCY);

	//set up the table [W^0, W^1, W^2, ..., W^(n/2 - 1)] where W = e^(-2*pi*i/n) for forward fft and W = e^(2*pi*i/n) for inverse fft
	double complex *exptable = malloc((SIZE / 2) * sizeof(double complex));
	if (exptable == NULL)
		return false;
	for (size_t i = 0; i < SIZE / 2; i++)
		exptable[i] = cexp(2 * M_PI * i / SIZE * I);

	DFT(in, out, SIZE);

    for(int i = 0; i < SIZE; ++i){
        printf("%.2lf j%.2lf   ", creal(in[i]), cimag(in[i]));
    }
    printf("\n");
    printf("\n");
    for(int i = 0; i < SIZE; ++i){
        printf("%.2lf j%.2lf   ", creal(out[i]), cimag(out[i]));
    }
    printf("\n");

    bool ret = fft_r2(inp, SIZE, exptable);
    printf("\n");
    for(int i = 0; i < SIZE; ++i){
        printf("%.2lf j%.2lf   ", creal(inp[i]), cimag(inp[i]));
    }
    printf("\n");

	free(exptable);
	return 0;
}