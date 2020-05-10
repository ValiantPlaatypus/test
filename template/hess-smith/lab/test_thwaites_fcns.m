
clear all; clc; close all

m = -0.1:0.005:0.1;
lam = -m;

ell = thwaites_ell( lam );
H   = thwaites_H(   lam );

figure
plot(m, ell, 'r-', 'LineWidth', 2), grid on, xlabel('m'), ylabel('l')

figure
plot(m, H  , 'r-', 'LineWidth', 2), grid on, xlabel('m'), ylabel('H')
