clear; close; clc;
figure; hold on;
ylim([0, 20])

n = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 30];
res_1 = [7.73, 5.721, 5.178, 4.946, 4.821, 4.744, 4.694, 4.659, 4.635, 4.6618, 4.6021, 4.591151, 4.5824034, 4.5753101, 4.5695062, 4.5542622, 4.5390622];
plot(n, res_1, 'bx--');
res_2 = [NaN, 15.46, 8.241, 6.565, 5.859, 5.478, 5.241, 5.085, 4.976, 4.89682, 4.8368, 4.7903634, 4.7538632, 4.724505, 4.7004498, 4.6378876, 4.576247];
plot(n, res_2, 'rx--');
res_3 = [NaN, NaN, 29.564, 11.882, 8.457, 7.086, 6.358, 5.912, 5.617, 5.41094, 5.26029, 5.1455441, 5.0561286, 4.9853691, 4.9281933, 4.7810472, 4.6389304];
plot(n, res_3, 'gx--');
res_4 = [NaN, NaN, NaN, 57.318, 17.123, 10.92, 8.612, 7.43, 6.722, 6.25806, 5.93605, 5.7003848, 5.5207693, 5.3805543, 5.2692477, 4.9898905, 4.7282449];
plot(n, res_4, 'mx--');
res_5 = [NaN, NaN, NaN, NaN, 115.262, 24.899, 14.125, 10.476, 8.698, 7.66336, 6.99861, 6.5402192, 6.2058174, 5.9516487, 5.7535587, 5.274786, 4.8459199];
plot(n, res_5, 'cx--');

xlabel("Number of atoms in the lattice");
ylabel("Wavelength of the resonant wave");
title("Wavelength of a resonant wave incident on a 1d lattice of regularly spaced atomes (d=1)");
legend("resonnance 1", "resonnance 2", "resonnance 3", "resonnance 4", "resonnance 5")
copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 12);
set(gca, "fontsize", 12)