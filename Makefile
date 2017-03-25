main:
	gcc -fopenmp -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas Main.c -O2 -std=c99
clean:
	rm a.out
	rm Knot.e*
	rm Knot.o*

run:
	./a.out 8_19Params.txt

trace:
	instruments -t "Time Profiler" -D profile_results ./a.out 8_19Params.txt
