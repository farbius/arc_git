%%
%
clc
clear
close all

%% constants
m2sec = 1000/3600;
c     = 299792458;

LOAD_MAT = 0;

%% signal generator
Fs      = 500e6;        % ADC rate
Tm      = 10.0e-6;    % period
F2      = 50e6;        % start frequency
F1      = 0e6;     % stop frequency 
Ns      = ceil(Tm*Fs)
ts      = (0:Ns-1)./Fs;
devF    = F2 - F1;

Ta      = 0.01;
Na      = ceil(Ta/Tm)

ta      = (0:Na-1).*Tm;

s       = exp(2*1i*pi*(-devF/2.*ts + devF/Tm/2.*ts.^2));
    
pspectrum(s,Fs,'spectrogram','FrequencyLimits',[0 Fs],'TimeResolution',.5e-6);

%% targets parameters
xt1   = 500;
yt1   = 50;
vt1   = -90*m2sec;

xt2   = 400;
yt2   = 50;
vt2   = 0*m2sec;

%%
Rt1   = sqrt((xt1 - vt1.*ta).^2 + yt1^2);
Rt2   = sqrt((xt2 - vt2.*ta).^2 + yt2^2);
s_raw = zeros(Na, Ns);

if LOAD_MAT == 1
    load('raw.mat');
else
    for i = 1 : Na

        td1          = 2*Rt1(i)/c; % time delay every Tm modulation period
        td2          = 2*Rt2(i)/c; % time delay every Tm modulation period
        s_raw(i, :) = exp(-2*1i*pi*(-devF/2.*(ts-td1) + devF/Tm/2.*(ts-td1).^2)).*exp(1i*4*pi*Rt1(i)/0.03) + ...
            + exp(-2*1i*pi*(-devF/2.*(ts-td2) + devF/Tm/2.*(ts-td2).^2)).*exp(1i*4*pi*Rt2(i)/0.03);  
    end
    
    save('raw.mat', 's_raw');
end

figure
plot(ts./1e-6, real(s_raw(1, :)))
title('period of raw signal: time domain')
xlabel('ts, usec')
grid on

figure
imagesc(real(s_raw))
title('raw signal: time domain')
xlabel('ts, usec')
ylabel('ta, sec')
grid on

%% demodulation Range dir
xop   = 400;
yop   = 50;
Rop   = sqrt(xop.^2 + yop^2);
top   = 2*Rop/c; 
s_op  = exp(2*1i*pi*(-devF/2.*(ts - top) + devF/Tm/2.*(ts - top).^2));
s_rng = zeros(Na, Ns);
for i = 1 : Na
    s_rng(i, :) = ifft(fft(s_raw(i, :)).*(fft(s_op)));
end

figure
imagesc(abs(s_rng))
title('FM demodulation')
xlabel('Range samples')
ylabel('Dopler samples')
grid on

figure
plot(abs(s_rng(1, :)))
grid on

figure
plot(real(s_rng(:, 1344)))
grid on

figure
plot(abs(fftshift(fft(s_rng(:, 1344)))))
grid on

%% dopler FFT
s_dpl = zeros(Na, Ns);
for i = 1 : Ns
    s_dpl(:, i) = fftshift(fft(s_rng(:, i)));
end

figure
imagesc(abs(s_dpl))
title('Dopler FFT')
xlabel('Range samples')
ylabel('Dopler samples')
grid on


  
  
  
