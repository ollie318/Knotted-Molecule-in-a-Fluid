main:
	gcc Main.c -O2

clean:
	rm a.out

run:
	./a.out 8_19Params.txt

trace:
	instruments -t "Time Profiler" -D profile_results ./a.out 8_19Params.txt
