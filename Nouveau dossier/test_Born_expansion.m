clear; close; clc;

k = 1;
x1 = -1.64;
y1 = 2.64;

theta = 0:0.01:2*pi;
r = 1;
x = r*cos(theta);
y = r*sin(theta);

figure; hold on;
A = -1i/sqrt(8*pi*k);
res = A * exp(i*k*r)/sqrt(r);
for n = 1:4
  born = A^n / ((x1^2 + y1^2)^(n/4));
  res = res +A * exp(i*k*r)/sqrt(r) * exp(-1i*k*(cos(theta)*x1 + sin(theta)*y1) + 1i*k*sqrt(x1^2 + y1^2)) * born;
  plot(theta, abs(res).^2);
endfor
legend("order=1", "order=2", "order=3");
xlabel("\\theta [rad]", 'interpreter', 'tex')
ylabel("|\\phi(\\theta)|²", 'interpreter', 'tex')
%ylim([0, 0.02]);
xlim([0, 2*pi]);
title("Born expansion of the wave function scattered by a delta potential")