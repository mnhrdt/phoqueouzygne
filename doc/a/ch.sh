#!/bin/sh

cat <<END | octave
K = 2;
T = 5;
N = 1001;
x = linspace(-3*T, 3*T, N);
dx = (x(2)-x(1));
h = exp(i*K*x.*x) .* (abs(x)<T);
f = conv(h, conj(h), "same") * dx;
s = (abs(x)<2*T) .* (2*T-abs(x)) .* sinc(K*x.*(2*T-abs(x))/pi);
xhfs = [x ; real(h) ; imag(h) ; real(f) ; imag(f) ; s ];
dlmwrite("o/xhfs", xhfs', " ");
END

cat <<END | gnuplot -persist
set zeroaxis
set title "N=101"
sinc(x)=sin(x)/x
K=2
T=5
plot [-12:12] [-2.5:11] \
	"o/xhfs" using 1:4 w points title "Re(h*h)" lt 2, \
	((abs(x)<2*T)*(2*T-abs(x))*sinc(K*x*(2*T-abs(x))))**1 w lines title "s" lt 1
END

#$paste o/x.txt o/hr.txt > o/hr
#$paste o/x.txt o/hi.txt > o/hi
#$paste o/x.txt o/fr.txt > o/fr
#$paste o/x.txt o/fi.txt > o/fi
#$paste o/x.txt o/s.txt > o/s
#
#cat <<END | gnuplot -persist
#plot [:] [-3:12] "o/hr" w lines, "o/hi" w lines, "o/fr" w lines, "o/fi" w lines, "o/s" w points pt 1
#END

#head o/xhfs
