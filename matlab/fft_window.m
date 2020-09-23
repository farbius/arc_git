clc
clear
close all


Fs  = 200e6;
f0  = 10e6;
tau = 50e-6;
Tp  = 1e-3;
Np  = 100;

N   = tau*Fs;
Ntp =  Tp*Fs
t   = [0:N-1]./Fs;
s   = cos(2*pi*f0.*t);

sp  = [s zeros(1,Ntp-N)];

x   = repmat(sp, [1 Np]);
x   = [zeros(1, ceil(1000*rand(1))) x];

% figure
% plot(x, '.-b')
% ylim([-1.5 1.5])
% grid on

Nwindow = 8192;

Nstep = ceil(length(x)/Nwindow) - 1;

Yf    = zeros(Nwindow/2, Nstep);

for k = 0:Nstep-1
    
    samples    = x(Nwindow*k+1:Nwindow*(k+1));
    Sf         = abs(fft(samples));
    Yf(:, k+1) = Sf(1:Nwindow/2);   
end


figure
imagesc(Yf)
grid on

figure
plot(Yf(411, :))
grid on

smpl = Yf(411, :);
I = zeros(1, Nstep);
I(smpl > 2000) = 1;

II = find(I > 0);
III = diff(II);






