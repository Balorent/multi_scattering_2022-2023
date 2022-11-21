clear; close; clc;

figure;
module = 0:0.001:1;
f  = (1.0 - 1.0 ./ (1 + module.^2)) .^ 0.2;
plot(module, f);

figure;
module = 0:0.001:1;
g = 0.1;
ctr = 0.1;
f  = 1/pi * atan2(g/2, (ctr - module));
plot(module, f);