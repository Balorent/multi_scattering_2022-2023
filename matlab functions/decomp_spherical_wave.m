clear; close; clc;
k0 = 2*pi + 0.001*1i;

figure;
[x, y] = meshgrid(-10:.1:10);
f = @(k) 1i/(4*pi) * exp(1i*k*x + 1i*sqrt(k0^2 - k^2)*abs(y)) / sqrt(k0^2 - k^2);
res = integral(f, -10, 10, 'ArrayValued', true);
mesh(x, y, abs(res));

figure;
r = 0.01:0.01:10;
f2 = @(k) 1/(1i*4*pi) * exp(1i*k*r) / sqrt(k0^2 - k^2);
res2 = integral(f2,  -20, 20,  'ArrayValued',  true);
plot(r, arg(res2));
hold on;
green = @(k, r) -1/(2*pi) * besselk(0, -1i*k*r, 0);
plot(r, arg(green(k0, r)));