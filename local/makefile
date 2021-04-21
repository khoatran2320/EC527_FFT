fft: serial_fft.c
	gcc -std=c99 serial_fft.c -lm -o fft 

fft2d: fft2d.c
	gcc -std=c99 fft2d.c -lm -o fft2d 

mtfft: mt_fft.c
	gcc -std=c99 mt_fft.c -lm -lpthread -o mtfft
.phony: clean

clean:
	rm -f fft mtfft fft2d