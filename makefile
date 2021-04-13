fft: serial_fft.c
	gcc -std=c99 serial_fft.c -lm -o fft

.phony: clean

clean:
	rm -f fft