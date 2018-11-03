#!/bin/sh

cat <<END | octave
K = 0.6;
T = 5;
N = 101;
x = linspace(-3*T, 3*T, N);
dx = (x(2)-x(1));
h = exp(i*K*x.*x) .* (abs(x)<T);
f = conv(h, conj(h), "same") * dx;
s = (abs(x)<2*T) .* sinc(K*x.*(2*T-abs(x))/pi);
s = (abs(x)<2*T) .* (2*T-abs(x)) .* sinc(K*x.*(2*T-abs(x))/pi);
dlmwrite("o/x.txt", x, "\n");
dlmwrite("o/hr.txt", real(h), "\n");
dlmwrite("o/hi.txt", imag(h), "\n");
dlmwrite("o/fr.txt", real(f), "\n");
dlmwrite("o/fi.txt", imag(f), "\n");
dlmwrite("o/s.txt", s, "\n");
END

paste o/x.txt o/hr.txt > o/hr
paste o/x.txt o/hi.txt > o/hi
paste o/x.txt o/fr.txt > o/fr
paste o/x.txt o/fi.txt > o/fi
paste o/x.txt o/s.txt > o/s

cat <<END | gnuplot -persist
plot [:] [-3:12] "o/hr" w lines, "o/hi" w lines, "o/fr" w lines, "o/fi" w lines, "o/s" w points pt 1
END
