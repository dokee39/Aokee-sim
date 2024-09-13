%% forward_solve
l_a = 60;
l_b = 150;
l_c = 220;

varphi_1 = pi * 2 / 3;
varphi_2 = pi * 3 / 5;
% [varphi_1, varphi_2] = meshgrid((pi / 2):0.05:(pi * 4 / 3), (pi * 4 / 3):-0.05:(pi / 2));

x_B_1 = l_a - l_b .* cos(varphi_1);
x_B_2 = -l_a + l_b .* cos(varphi_2);
y_B_1 = l_b .* sin(varphi_1);
y_B_2 = l_b .* sin(varphi_2);

x_B = x_B_1 - x_B_2;
y_B = y_B_1 - y_B_2;
varphi_B_1 = acos(sqrt(x_B .* x_B + y_B .* y_B) / (2 .* l_c)) - atan2(y_B, x_B);

x_C = l_a - l_b .* cos(varphi_1) - l_c .* cos(varphi_B_1); 
y_C = l_b .* sin(varphi_1) + l_c .* sin(varphi_B_1);

% figure;
% scatter(x_C, y_C);

l = sqrt(x_C .* x_C + y_C .* y_C);
theta_l = atan2(x_C, y_C);

% figure;
% surf(varphi_1, varphi_2, l);
% figure;
% surf(varphi_1, varphi_2, theta_l);

%% inverse_solve
l_1 = sqrt(l_a .* l_a + l .* l - 2 .* l_a .* l .* sin(theta_l)); 
l_2 = sqrt(l_a .* l_a + l .* l + 2 .* l_a .* l .* sin(theta_l)); 

varphi_1_new = acos((l_1 .* l_1 + l_a .* l_a - l .* l) / (2 .* l_1 .* l_a)) + acos((l_1 .* l_1 + l_b .* l_b - l_c .* l_c) / (2 .* l_1 .* l_b));
varphi_2_new = acos((l_2 .* l_2 + l_a .* l_a - l .* l) / (2 .* l_2 .* l_a)) + acos((l_2 .* l_2 + l_b .* l_b - l_c .* l_c) / (2 .* l_2 .* l_b));

%% jacobian matrix
% SJTU
one_over_l_b = 1 ./ l_b;
l_over_l_b = l ./ l_b;

x_CB_1 = x_C - x_B_1;
y_CB_1 = y_C - y_B_1;
x_CB_2 = x_C - x_B_2;
y_CB_2 = y_C - y_B_2;

sin_theta_l = sin(theta_l);
cos_theta_l = cos(theta_l);
sin_varphi_1 = sin(varphi_1);
cos_varphi_1 = cos(varphi_1);
sin_varphi_2 = sin(varphi_2);
cos_varphi_2 = cos(varphi_2);

x_CB_1_times_sin_varphi_1 = x_CB_1 .* sin_varphi_1;
y_CB_1_times_cos_varphi_1 = y_CB_1 .* cos_varphi_1;
x_CB_2_times_sin_varphi_2 = x_CB_2 .* sin_varphi_2;
y_CB_2_times_cos_varphi_2 = y_CB_2 .* cos_varphi_2;

J_inverse = [(one_over_l_b .* (x_CB_1 .* sin(theta_l) + y_CB_1 .* cos(theta_l)) ./ (x_CB_1_times_sin_varphi_1 + y_CB_1_times_cos_varphi_1)) ...
             (l_over_l_b .* (x_CB_1 .* cos(theta_l) - y_CB_1 .* sin(theta_l)) ./ (x_CB_1_times_sin_varphi_1 + y_CB_1_times_cos_varphi_1)); ...
             (- one_over_l_b .* (x_CB_2 .* sin(theta_l) + y_CB_2 .* cos(theta_l)) ./ (x_CB_2_times_sin_varphi_2 - y_CB_2_times_cos_varphi_2)) ...
             (- l_over_l_b .* (x_CB_2 .* cos(theta_l) - y_CB_2 .* sin(theta_l)) ./ (x_CB_2_times_sin_varphi_2 - y_CB_2_times_cos_varphi_2))];
J = inv(J_inverse)

% HEU
l_b_over_l = l_b / l;
varphi_B_2 = atan2(y_C - y_B_2, x_C - x_B_2);
sin_varphi_B_1_plus_varphi_B_2 = sin(varphi_B_1 + varphi_B_2);

J = [(-l_b .* cos(theta_l + varphi_B_2) .* sin(varphi_1 - varphi_B_1) ./ sin_varphi_B_1_plus_varphi_B_2) ...
     (-l_b .* cos(theta_l - varphi_B_1) .* sin(varphi_2 - varphi_B_2) ./ sin_varphi_B_1_plus_varphi_B_2); ...
     (l_b_over_l .* sin(theta_l + varphi_B_2) .* sin(varphi_1 - varphi_B_1) ./ sin_varphi_B_1_plus_varphi_B_2) ...
     (l_b_over_l .* sin(theta_l - varphi_B_1) .* sin(varphi_2 - varphi_B_2) ./ sin_varphi_B_1_plus_varphi_B_2)]

