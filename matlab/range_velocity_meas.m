clc
clear

c    = 299792458;
f0   = 2.3e9;
la   = c/f0;

fdev = 20e6;
TCPI = 100e-3;


% to be measured in real world application
fB1 = 2000;
fB2 = -500;

% matrix constant
a = 2*fdev/(c*TCPI);
b = 2/la;

% coeff matrix
A = [a -b; -a -b];

B = [fB1; fB2];

% solve linear equations:
x = A\B;

R  = x(1)
vr = x(2)




