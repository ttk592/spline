# adjust path to ALGLIB library
LALGLIB=-L../../alglib/ -lalglib
IALGLIB=-I../../alglib/src

# compiler setting
CC=g++
CFLAGS = ${IALGLIB} -I../src
CFLAGS += -O2 -DNDEBUG
CFLAGS += -Wall
LFLAGS=

PROGS=demo interpol interpol_alglib bench

all:		${PROGS}

demo:		demo.o
		${CC} $< -o $@

interpol:	interpol.o
		${CC} $< -o $@

interpol_alglib:	interpol_alglib.o
		${CC} $< -o $@ ${LALGLIB}

bench:		bench.o
		${CC} $< -o $@ ${LALGLIB}

plot:		interpol interpol_alglib
		./interpol > out.txt
		./interpol_alglib >> out.txt
		gnuplot plot.txt

%.o: %.cpp
		${CC} ${CFLAGS} -c $< -o $@

clean:
		rm -f *.o ${PROGS} out.txt
