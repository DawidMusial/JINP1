aprox: main.o splines.o points.o  aproksymator_f_trygonometrycznej.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o  aproksymator_f_trygonometrycznej.o -L gaus -lm -l ge -Wall -ansi -pedantic 

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

aproksymator_f_trygonometrycznej.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_f_trygonometrycznej.c 

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta
