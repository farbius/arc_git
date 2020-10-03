clear
clc
close all
N = 100;  % number of samples
a = 0.1; % acceleration
sigmaPsi=1;
sigmaEta=50;
k=1:N;
x=k;
x(1)=0;
z(1)=x(1)+normrnd(0,sigmaEta);

for t=1:(N-1)
   x(t+1)= x(t)+a*t+normrnd(0,sigmaPsi); 
   z(t+1)= x(t+1)+normrnd(0,sigmaEta);
end

x = repmat(x, [1 2]);
z = repmat(z, [1 2]);
k=1:2*N;
%kalman filter
xOpt(1)=z(1);
eOpt(1)=sigmaEta; % eOpt(t) is a square root of the error dispersion (variance). It's not a random variable. 
for t=1:(2*N-1)
  eOpt(t+1)=sqrt((sigmaEta^2)*(eOpt(t)^2+sigmaPsi^2)/(sigmaEta^2+eOpt(t)^2+sigmaPsi^2));
  K(t+1)=(eOpt(t+1))^2/sigmaEta^2;
 xOpt(t+1)=(xOpt(t)+a*t)*(1-K(t+1))+K(t+1)*z(t+1);
end
plot(k,xOpt,'.-r',k,z,'.-b',k,x,'.-')
legend('xOpt','z','x_real')
grid on