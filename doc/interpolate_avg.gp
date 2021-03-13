set term pdfcairo color solid size  12cm,8cm
set output "interpolate_avg.pdf"

# input points to interpolate (0,y1), (h,y2) with average avg
h=1.0
avg=1.0
y1=4.0
y2=2.0

# calculate quadratic function: f(x) = a + b*x + c*x^2
# which has the average avg and goes through both points
a=y1
b=2.0/h*(-2.0*y1-y2+3.0*avg)
c=3.0/(h*h)*(y1+y2-2.0*avg)
f(x) = a + b*x + c*x*x

# make positive
R=sqrt( (y1/avg)**2 + (y2/avg)**2 )
if(R>3) {
    z1=y1/R*3.0
    z2=y2/R*3.0
} else {
    z1=y1
    z2=y2
}
aa=z1
bb=2.0/h*(-2.0*z1-z2+3.0*avg)
cc=3.0/(h*h)*(z1+z2-2.0*avg)
g(x) = aa + bb*x + cc*x*x

#set size ratio -1
set samples 300
set colors classic
set xzeroaxis
set yzeroaxis
set xrange [-0.05*h:1.05*h]

set arrow from h, graph 0 to h,graph 1 nohead lt 0

set arrow from 0,y1 to 0,z1 head lt 2 lw 1 
set arrow from h,y2 to h,z2 head lt 2 lw 1
set label " =h" at h,screen 0.05
set label "  y_1" at 0,y1
set label "  y_2" at h,y2


plot	f(x) notitle with l lt 1 lw 1,\
	g(x) notitle with l lt 2 lw 1 dt 2,\
	avg  notitle with l lt 0 lw 1,\
	"-" using 1:(f($1)) notitle with p pt 7 lt 1,\
	"-" using 1:(g($1)) notitle with p pt 7 lt 2
0.0
1.0 
e
0.0
1.0 
e


pause -1
