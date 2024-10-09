%% 
l_a = 0.075;
l_b = 0.15;
l_c = 0.22:0.001:0.35;

varphi_1 = pi * 2 / 3;
varphi_2 = pi * 2 / 3;

x_B_1 = l_a - l_b .* cos(varphi_1);
x_B_2 = -l_a + l_b .* cos(varphi_2);
y_B_1 = l_b .* sin(varphi_1);
y_B_2 = l_b .* sin(varphi_2);

x_B = x_B_1 - x_B_2;
y_B = y_B_1 - y_B_2;
varphi_B_1 = acos(sqrt(x_B .* x_B + y_B .* y_B) ./ (2 .* l_c)) - atan2(y_B, x_B);

x_C = l_a - l_b .* cos(varphi_1) - l_c .* cos(varphi_B_1); 
y_C = l_b .* sin(varphi_1) + l_c .* sin(varphi_B_1);

l = sqrt(x_C .* x_C + y_C .* y_C);
theta = atan2(x_C, y_C);

l_b_over_l = l_b ./ l;
varphi_B_2 = atan2(y_C - y_B_2, x_C - x_B_2);
sin_varphi_B_1_plus_varphi_B_2 = sin(varphi_B_1 + varphi_B_2);

f_1 = (-l_b .* cos(theta + varphi_B_2) .* sin(varphi_1 - varphi_B_1) ./ sin_varphi_B_1_plus_varphi_B_2);
f_2 = (-l_b .* cos(theta - varphi_B_1) .* sin(varphi_2 - varphi_B_2) ./ sin_varphi_B_1_plus_varphi_B_2);
tau_1 = (l_b_over_l .* sin(theta + varphi_B_2) .* sin(varphi_1 - varphi_B_1) ./ sin_varphi_B_1_plus_varphi_B_2);
tau_2 = (l_b_over_l .* sin(theta - varphi_B_1) .* sin(varphi_2 - varphi_B_2) ./ sin_varphi_B_1_plus_varphi_B_2);

figure;
plot(l_c, f_1, l_c, f_2, l_c, tau_1, l_c, tau_2);
title("l_c");
legend("f_1", "f_2", "tau_1", "tau_2");

%% 
l_a = 0.075;
l_b = 0.10:0.001:0.20;
l_c = 0.27;

varphi_1 = pi * 2 / 3;
varphi_2 = pi * 2 / 3;

x_B_1 = l_a - l_b .* cos(varphi_1);
x_B_2 = -l_a + l_b .* cos(varphi_2);
y_B_1 = l_b .* sin(varphi_1);
y_B_2 = l_b .* sin(varphi_2);

x_B = x_B_1 - x_B_2;
y_B = y_B_1 - y_B_2;
varphi_B_1 = acos(sqrt(x_B .* x_B + y_B .* y_B) ./ (2 .* l_c)) - atan2(y_B, x_B);

x_C = l_a - l_b .* cos(varphi_1) - l_c .* cos(varphi_B_1); 
y_C = l_b .* sin(varphi_1) + l_c .* sin(varphi_B_1);

l = sqrt(x_C .* x_C + y_C .* y_C);
theta = atan2(x_C, y_C);

l_b_over_l = l_b ./ l;
varphi_B_2 = atan2(y_C - y_B_2, x_C - x_B_2);
sin_varphi_B_1_plus_varphi_B_2 = sin(varphi_B_1 + varphi_B_2);

f_1 = (-l_b .* cos(theta + varphi_B_2) .* sin(varphi_1 - varphi_B_1) ./ sin_varphi_B_1_plus_varphi_B_2);
f_2 = (-l_b .* cos(theta - varphi_B_1) .* sin(varphi_2 - varphi_B_2) ./ sin_varphi_B_1_plus_varphi_B_2);
tau_1 = (l_b_over_l .* sin(theta + varphi_B_2) .* sin(varphi_1 - varphi_B_1) ./ sin_varphi_B_1_plus_varphi_B_2);
tau_2 = (l_b_over_l .* sin(theta - varphi_B_1) .* sin(varphi_2 - varphi_B_2) ./ sin_varphi_B_1_plus_varphi_B_2);

figure;
plot(l_b, f_1, l_b, f_2, l_b, tau_1, l_b, tau_2);
title("l_b");
legend("f_1", "f_2", "tau_1", "tau_2");

%% 
l_a = 0.03:0.001:0.10;
l_b = 0.15;
l_c = 0.27;

varphi_1 = pi * 2 / 3;
varphi_2 = pi * 2 / 3;

x_B_1 = l_a - l_b .* cos(varphi_1);
x_B_2 = -l_a + l_b .* cos(varphi_2);
y_B_1 = l_b .* sin(varphi_1);
y_B_2 = l_b .* sin(varphi_2);

x_B = x_B_1 - x_B_2;
y_B = y_B_1 - y_B_2;
varphi_B_1 = acos(sqrt(x_B .* x_B + y_B .* y_B) ./ (2 .* l_c)) - atan2(y_B, x_B);

x_C = l_a - l_b .* cos(varphi_1) - l_c .* cos(varphi_B_1); 
y_C = l_b .* sin(varphi_1) + l_c .* sin(varphi_B_1);

l = sqrt(x_C .* x_C + y_C .* y_C);
theta = atan2(x_C, y_C);

l_b_over_l = l_b ./ l;
varphi_B_2 = atan2(y_C - y_B_2, x_C - x_B_2);
sin_varphi_B_1_plus_varphi_B_2 = sin(varphi_B_1 + varphi_B_2);

f_1 = (-l_b .* cos(theta + varphi_B_2) .* sin(varphi_1 - varphi_B_1) ./ sin_varphi_B_1_plus_varphi_B_2);
f_2 = (-l_b .* cos(theta - varphi_B_1) .* sin(varphi_2 - varphi_B_2) ./ sin_varphi_B_1_plus_varphi_B_2);
tau_1 = (l_b_over_l .* sin(theta + varphi_B_2) .* sin(varphi_1 - varphi_B_1) ./ sin_varphi_B_1_plus_varphi_B_2);
tau_2 = (l_b_over_l .* sin(theta - varphi_B_1) .* sin(varphi_2 - varphi_B_2) ./ sin_varphi_B_1_plus_varphi_B_2);

figure;
plot(l_a, f_1, l_a, f_2, l_a, tau_1, l_a, tau_2);
title("l_a");
legend("f_1", "f_2", "tau_1", "tau_2");

