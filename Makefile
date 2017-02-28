main:
	gcc -fprofile-arcs -ftest-coverage Main.c -O2 -pg

clean:
	rm a.out

trace:
	instruments -t "Time Profiler" -D profile_results ./a.out 8_19Params.txt
