clear; close; clc;

G = @(k, r) -1/(2*pi) * besselk (0, -i*k*r, 0);
G2 = @(k, r) exp(i*k*r)./r;
G3 = @(k, r) -sqrt(-1./(8*pi*i*k*r)) .* exp(i*k*r);
k = 2*pi;

figure; hold on;
r = 0.001:0.0001:5;
plot(r, arg(G(k, r)), 'b');
plot(r, arg(G2(k, r)), 'r');
legend("green (k=2*pi)", "exp(ikr)/r (k=2*pi)");
title("Comparaison des phases")

figure; hold on;
r = 0.001:0.0001:10;
plot(r, abs(G(k, r)), 'b');
plot(r, abs(G2(k, r)), 'r');
legend("green (k=2*pi)", "exp(ikr)/r (k=2*pi)");
title("Comparaison des modules")
ylim([0, 1])

figure; hold on;
r = 0.001:0.0001:2;
plot(r, arg(G(k, r)), 'b');
plot(r, arg(G3(k, r)), 'r')
legend("green (k=2*pi)", "approx (k=2*pi)");
title("Comparaison des phases")

figure; hold on;
r = 0.001:0.0001:1;
plot(r, abs(G(k, r)), 'b');
plot(r, abs(G3(k, r)), 'r');
legend("green (k=2*pi)", "approx (k=2*pi)");
title("Comparaison des modules")
ylim([0, 1]);