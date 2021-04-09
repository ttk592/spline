# some programmes require ALGLIB for comparison
# adjust the path to the ALGLIB library below
ALGLIB=../../alglib
LALGLIB=-L$(ALGLIB) -lalglib
IALGLIB=-I$(ALGLIB)/src

# compiler setting
CC=g++
CFLAGS = ${IALGLIB} -I../src -std=c++11 -Wall -Wextra
CFLAGS += -O2
#CFLAGS += -DNDEBUG			# disables all assert() checks
CFLAGS += -DHAVE_SSTREAM		# whether compiler has <sstream>
#CFLAGS +=-D_GLIBCXX_USE_CXX11_ABI=0	# for new gcc & boost compiled with old
LFLAGS=


PROGS=simple_demo plot plot_avg_preserv plot_parametric plot_bspline cubic_root
PROGS_ALGLIB=plot_alglib bench
PROGS_TESTS=tests/unit_tests tests/compile_units tests/simple_demo_c98

all:		${PROGS}
more:		${PROGS_ALGLIB}
tests:		${PROGS_TESTS}
		./tests/unit_tests


simple_demo:	simple_demo.o
		${CC} $< -o $@
plot:		plot.o
		${CC} $< -o $@
plot_avg_preserv:	plot_avg_preserv.o
		${CC} $< -o $@
plot_parametric:	plot_parametric.o
		${CC} $< -o $@
plot_bspline:	plot_bspline.o
		${CC} $< -o $@
cubic_root:	cubic_root.o gsl_solve_cubic.o
		${CC} $^ -o $@

tests/compile_units:	tests/compile_units.o tests/myfunc1.o tests/myfunc2.o
		${CC} $? -o $@ 
tests/simple_demo_c98:	tests/simple_demo_c98.cpp ../src/spline.h
		${CC} -Wall -Wextra -Werror -I../src -std=c++98 $< -o $@
tests/unit_tests:	tests/unit_tests.o
		${CC} $< -o $@ -lboost_unit_test_framework


%.o: 		%.cpp ../src/spline.h
		${CC} ${CFLAGS} -c $< -o $@

tests/%.o:	tests/%.cpp ../src/spline.h
		${CC} ${CFLAGS} -Werror -c $< -o $@


check_alglib:
ifeq ($(wildcard $(ALGLIB)),)
	$(error Please install ALGLIB and adjust the corresponding variables within this Makefile)
endif

plot_alglib:	check_alglib plot_alglib.o
		${CC} plot_alglib.o -o $@ ${LALGLIB}
bench:		check_alglib bench.o
		${CC} bench.o -o $@ ${LALGLIB}


clean:
		rm -f *.o tests/*.o ${PROGS} ${PROGS_ALGLIB} ${PROGS_TESTS}  *.csv
