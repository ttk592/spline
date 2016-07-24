# some programmes require ALGLIB for comparison
# adjust the path to the ALGLIB library below
ALGLIB=../../alglib
LALGLIB=-L$(ALGLIB) -lalglib
IALGLIB=-I$(ALGLIB)/src

# compiler setting
CC=g++
CFLAGS = ${IALGLIB} -I../src
CFLAGS += -O2 -DNDEBUG
CFLAGS += -Wall
LFLAGS=



PROGS=demo interpol check_alglib interpol_alglib bench

all:		${PROGS}

demo:		demo.o
		${CC} $< -o $@

interpol:	interpol.o
		${CC} $< -o $@

check_alglib:
ifeq ($(wildcard $(ALGLIB)),)
  $(error Please install ALGLIB and adjust the corresponding variables within this Makefile)
endif

interpol_alglib:	interpol_alglib.o
		${CC} $< -o $@ ${LALGLIB}

bench:		bench.o
		${CC} $< -o $@ ${LALGLIB}


%.o: %.cpp
		${CC} ${CFLAGS} -c $< -o $@

clean:
		rm -f *.o ${PROGS} *.csv
