clc; clear all; close all;

% 1.a.)

c = 1;
t = 1:5;
x = -100:.1:100;
%u_xt = @(x, t) exp(-(x - c * t).^2); 
 
%u_xt2 = @(x, t) exp(-(x - c * t^2).^2); % part c.)

u_xt = @(x, t) (1/t) * exp(-(x - c * t).^2);  % part d.)

f_list = zeros(length(t), length(x));
f_list2 = zeros(length(t), length(x));

for n = 1:length(t)
    f_list(n, :) = u_xt(x, t(n));


end

%%

for n = 1:5
   plot(x, f_list(n, :))
   xlim([-5 50])
   ylim([-.2 1.2])
   pause(1) 
   
end

[M, I] = max(f_list(4, :));
[M2, I2] = max(f_list(5, :));
speed = x(I2) - x(I);
speed
