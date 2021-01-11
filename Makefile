CC = gcc
OBJ: MM_parallel MM_serial strassen_parallel random strassen_serial strassen_hybrid_serial strassen_hybrid_parallel_section strassen_hybrid_parallel_for

%: %.c
	$(CC) -O3 -g -o $@ $< -fopenmp

run: $(OBJ)
	@./MM_serial
	@./MM_parallel
	@./strassen_serial
	@./strassen_parallel
	@./strassen_hybrid_serial
	@./strassen_hybrid_parallel_section
	@./strassen_hybrid_parallel_for 

	@echo "difference between MM_serial.txt and MM_parallel.txt :"
	@diff MM_serial.txt MM_parallel.txt

	@echo "difference between serial.txt and strassen_serial.txt :"
	@diff MM_serial.txt strassen_serial.txt

	@echo "difference between serial.txt and strassen_parallel.txt :"
	@diff MM_serial.txt strassen_parallel.txt

	@echo "difference between serial.txt and strassen_hybrid_serial.txt :"
	@diff MM_serial.txt strassen_hybrid_serial.txt

	@echo "difference between serial.txt and strassen_hybrid_parallel_section.txt :"
	@diff MM_serial.txt strassen_hybrid_parallel_section.txt

	@echo "difference between serial.txt and strassen_hybrid_parallel_for.txt :"
	@diff MM_serial.txt strassen_hybrid_parallel_for.txt

.PHONY: clean
clean:
	-rm -f MM_parallel MM_serial strassen_parallel strassen_serial strassen_hybrid_serial strassen_hybrid_parallel_section strassen_hybrid_parallel_for random 
clean-data:
	-rm -f *.txt
