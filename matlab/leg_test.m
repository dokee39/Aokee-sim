l_a = 60;
l_b = 150;
l_c = 220;

[varphi_1, varphi_2] = meshgrid((pi / 3):0.05:(pi * 4 / 3), (pi * 4 / 3):-0.05:(pi/ 2));

x_B_1 = l_a - l_b .* cos(varphi_1);
x_B_2 = -l_a + l_b .* cos(varphi_2);
y_B_1 = l_b .* sin(varphi_1);
y_B_2 = l_b .* sin(varphi_2);

x_B = x_B_1 - x_B_2;
y_B = y_B_1 - y_B_2;
varphi_B_1 = acos(sqrt(x_B .* x_B + y_B .* y_B) / (2 .* l_c)) - atan2(y_B, x_B);

x_C = l_a - l_b .* cos(varphi_1) - l_c .* cos(varphi_B_1); 
y_C = l_b .* sin(varphi_1) + l_c .* sin(varphi_B_1);

l = sqrt(x_C .* x_C + y_C .* y_C);
theta_l = atan2(x_C, y_C);

figure;
surf(varphi_1, varphi_2, l);
figure;
surf(varphi_1, varphi_2, theta_l);

l_1 = sqrt(l_a .* l_a + l .* l - 2 .* l_a .* l .* sin(theta_l)); 
l_2 = sqrt(l_a .* l_a + l .* l + 2 .* l_a .* l .* sin(theta_l)); 

varphi_1_new = acos((l_1 .* l_1 + l_a .* l_a - l .* l) / (2 .* l_1 .* l_a)) + acos((l_1 .* l_1 + l_b .* l_b - l_c .* l_c) / (2 .* l_1 .* l_b));
varphi_2_new = acos((l_2 .* l_2 + l_a .* l_a - l .* l) / (2 .* l_2 .* l_a)) + acos((l_2 .* l_2 + l_b .* l_b - l_c .* l_c) / (2 .* l_2 .* l_b));


