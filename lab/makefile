fft: serial_fft.c
	gcc -std=c99 serial_fft.c -lm -o fftfft: serial_fft.c
	gcc -std=c99 serial_fft.c -lm -o fft 

fft2d: fft2d.c
	gcc -std=c99 fft2d.c -lm -o fft2d 

mtfft: mt_fft.c
	gcc -std=gnu99 mt_fft.c  -lm -pthread -o mtfft

mtfft2: mtfft2.c
	gcc -O3 -std=gnu99 mtfft2.c  -lm -pthread -o mtfft2

mtfft3: mtfft3.c
	gcc -std=gnu99 mtfft3.c -lm -pthread -o mtfft3

test: test.c
	gcc -std=gnu99 test.c -lm -pthread -o test
.phony: clean

clean:
	rm -f fft mtfft mtfft2 fft2d mtfft3 test