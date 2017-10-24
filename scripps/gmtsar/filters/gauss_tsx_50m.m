%
%   matlab script to design a gaussian filter
%   for tsx after 4 x 2 decimation so the 0.5 gain
%   is at 50 km wavelength
%
clear
sigx=1.5;
sigy=1.5;
nx=5;
ny=5;
x=((-nx/2:(nx/2-1))+.5)/sigx;
y=((-ny/2:(ny/2-1))+.5)/sigy;
x2=x.*x;
y2=(y.*y)';
allx=ones(ny,1)*x2;
ally=y2*ones(1,nx);
r2=allx+ally;
gauss=exp(-.5*r2);
save gauss_tsx_50m -ascii gauss
