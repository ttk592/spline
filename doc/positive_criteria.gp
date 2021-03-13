set term pdfcairo color solid size  8cm,8cm
set output "positive_criteria.pdf"


# implicit equation: x^2+y^2+xy-6(x+y)+9=0
f(x) = -0.5*(x-6.0) + sqrt( (-0.75*x+3.0)*x )
g(x) = -0.5*(x-6.0) - sqrt( (-0.75*x+3.0)*x )

# sufficient criteria: x^2+y^2=9
h(x) = sqrt(9.0-x*x)


set size ratio -1
set samples 500
set colors classic
set xlabel "y_1"
set ylabel "y_2"
set xrange [0:4]

set arrow from 0,0 to 0,3 nohead lt 1 lw 2

plot	f(x) title "exact" with l lt 1 lw 2,\
	g(x)*(x>=3) notitle with l lt 1 lw 2,\
	h(x)*(x<=3) title "sufficient" with l lt 2 dt 2


pause -1
