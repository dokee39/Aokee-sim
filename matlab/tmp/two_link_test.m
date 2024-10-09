%% param
l_1 = 0.19;
l_2 = 0.23;

% varphi_1 = pi * 1 / 4;
% varphi_2 = pi * 1 / 2;
[varphi_1, varphi_2] = meshgrid(0:0.05:(pi / 2), pi:-0.05:0);

%% forward_solve
x_B = l_1 .* sin(varphi_1) + l_2 .* sin(varphi_1 - varphi_2);
y_B = l_1 .* cos(varphi_1) + l_2 .* cos(varphi_1 - varphi_2);

figure;
scatter(x_B, y_B);

l = sqrt(x_B .* x_B + y_B .* y_B);
theta = atan2(x_B, y_B);

% figure;
% surf(varphi_1, varphi_2, l);
% figure;
% surf(varphi_1, varphi_2, theta);

%% inverse_solve
varphi_1_new = theta + acos((l .* l + l_1 .* l_1 - l_2 .* l_2) ./ (2 .* l .* l_1));
varphi_2_new = pi + acos((l_1 .* l_1 + l_2 .* l_2 - l .* l) ./ (2 .* l_1 .* l_2));

%% jacobian matrix
J = [l_1 * sin(theta - varphi_1) - l_2 .* sin(theta + varphi_1 - varphi_2) ...
     l_2 * sin(theta + varphi_1 - varphi_2); ...
     - l_2 ./ l .* cos(theta + varphi_1 - varphi_2) - l_1 ./ l .* cos(theta - varphi_1) ...
     l_2 ./ l .* cos(theta + varphi_1 - varphi_2)];

%% force test
theta = (-pi / 6):0.05:(pi / 6);
l = 0.28 ./ cos(theta);

varphi_1 = theta + acos((l .* l + l_1 .* l_1 - l_2 .* l_2) ./ (2 .* l .* l_1));
varphi_2 = pi + acos((l_1 .* l_1 + l_2 .* l_2 - l .* l) ./ (2 .* l_1 .* l_2));

f_1 = l_1 * sin(theta - varphi_1) - l_2 .* sin(theta + varphi_1 - varphi_2);
f_2 = l_2 * sin(theta + varphi_1 - varphi_2);
tau_1 = - l_2 ./ l .* cos(theta + varphi_1 - varphi_2) - l_1 ./ l .* cos(theta - varphi_1);
tau_2 = l_2 ./ l .* cos(theta + varphi_1 - varphi_2);

figure;
plot(theta, f_1, theta, f_2, theta, tau_1, theta, tau_2);
legend("f_1", "f_2", "tau_1", "tau_2");
