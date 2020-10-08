%% FM CW Radar detection and parameters measurement
% 08/10/20 A.Rostov aleksei.rostov@protonmail.com
clc
clear
close all


load('centr.mat');
fc = Centras;
 N = length(Centras);

% проверка на широкополосность
dF_bin     = abs(diff(fc));
treshold   = mean(dF_bin, 'omitnan');
ind        = find(dF_bin > treshold);

f0         = abs(fc(ind(2)) - fc(ind(1)+1))...
                    / (fc(ind(2)) - abs(fc(ind(2)) - fc(ind(1)+1))/2);



% widebn_ind = find(dF_bin < 2);


% if ~isempty(widebn_ind) && length(widebn_ind) > 0.8*N
% Сигнал широкополосный, если dF > f0/10
if(f0 < 0.1) 
    disp("узкополосный сигнал")
    return
end




DF     = sum(dF_bin(ind))/ length(ind)*100e3;
d_ind  = sum(abs(diff(ind)))/(length(ind)-1)*8192/800e6;

dev = DF/1e6;
Tp  = d_ind/1e-6;
F1  = fc(ind(1)+1)*100e3/1e6;
F2  = fc(ind(1))*100e3/1e6;


figure
plot(fc, '.-b')
grid on
ylabel('Fbins');
xlabel(sprintf('deviation=%3.2f MHz\n period=%4.2f us\n start freq=%3.2f MHz\n stop freq=%3.2f MHz \n',...
    dev,Tp,F1,F2),'Color','r');


disp('Application Done')