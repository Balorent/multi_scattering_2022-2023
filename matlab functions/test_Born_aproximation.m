clear; close; clc;

k = 10;
xn = [5];
yn = [5];
sn = [(1+1i)];

theta = 0:0.001:2*pi;
r = 1;
x = r*cos(theta);
y = r*sin(theta);

figure; hold on;
A = -1i/sqrt(8*pi*k);
res = A * exp(i*k*r)/sqrt(r);
for n = 1:size(xn)(2)
  res = res + A * exp(i*k*r)/sqrt(r) * A * sn(n) * exp(-1i*k*(cos(theta)*xn(n) + sin(theta)*yn(n)) + 1i*k*sqrt(xn(n)^2 + yn(n)^2)) / sqrt(xn(n)^2 + yn(n)^2);
endfor
plot(theta, abs(res).^2);
xlabel("\\theta [rad]", 'interpreter', 'tex')
ylabel("|\\phi(\\theta)|²", 'interpreter', 'tex')
%ylim([0, max(abs(res).^2)]);
xlim([0, 2*pi]);
title("Born expansion of the wave function scattered by a delta potential")