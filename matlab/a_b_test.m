clc
clear
close

N = 200;
dt = 1.0;
xk_1 = zeros(1, N); 
vk_1 = zeros(1, N); 
a = 0.85; 
b = 0.05;

% xm = randn(1, N);
load('cenrt.mat');
xm = Centras;

for i = 2 : length(xm)
    
    xk = xk_1(i-1) + (vk_1(i-1) * dt);
    vk = vk_1(i-1);
    
    rk = xm(i-1) - xk;
    
    xk = xk + a*rk;
    vk = (b * rk) / dt;
    
    vk_1(i) = vk;
    xk_1(i) = xk;
 
end

figure
plot(1:N, xm, '.-b', 1:N, xk_1, '.-r')
grid on



