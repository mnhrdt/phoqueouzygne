#!/bin/sh

cat <<END | octave
K = 1;
T = 5;
N = 1000;
x = linspace(-3*T, 3*T, N);
dx = (x(2)-x(1))*2;
h = cos(K*x.*x) .* (abs(x)<T);
f = conv(h, h, "same") * dx;
s = 2*(abs(x)<2*T) .* (5-abs(x)) .* sinc(K.*x.*(5-abs(x)/pi));
dlmwrite("o/x.txt", x, "\n");
dlmwrite("o/h.txt", h, "\n");
dlmwrite("o/f.txt", f, "\n");
dlmwrite("o/s.txt", s, "\n");
END

paste o/x.txt o/h.txt > o/h
paste o/x.txt o/f.txt > o/f
paste o/x.txt o/s.txt > o/s

cat <<END | gnuplot -persist
plot [-15:15] [-2:12] "o/h" w lines, "o/f" w lines, "o/s" w lines
END
