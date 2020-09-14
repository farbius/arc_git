function [xk_1, vk_1] = alpha_beta_filter(xm, xk_1, vk_1)
% Функция реализует стандартный альфа-бета-фильтр для
% предотвращения "забрасываний" порога ШП-помехами и сигналами
a_plus = 0.0001;
a_minus = 0.0008;
b = 0.000;
dt = 1;

xk = xk_1 + ( vk_1 * dt );
vk = vk_1;

rk = xm - xk;

if (rk > 0),
    xk = xk + a_plus * rk;
else
    xk = xk + a_minus*rk;
end

vk = vk + ( b * rk ) / dt;

xk_1 = xk;
vk_1 = vk;

