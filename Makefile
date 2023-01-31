all: a_star_algorithm

a_star_algorithm: a_star_algorithm.cpp
	mpic++ a_star_algorithm.cpp -o a_star_algorithm
.PHONY: clean
clean:
	rm -f $(ODIR)/*.o a_star_algorithm res.txt

run:
	mpirun -np 6 ./a_star_algorithm test/matrix.txt
show:
	cat res.txt
