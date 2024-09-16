syms l_a l_b l_c varphi_1(t) varphi_2(t) varphi_B_1(t) varphi_B_2(t) dot_varphi_1 dot_varphi_2 theta l

dot_x_B_1 = l_b * dot_varphi_1 * sin(varphi_1);
dot_x_B_2 = -l_b * dot_varphi_2 * sin(varphi_2);
dot_y_B_1 = l_b * dot_varphi_1 * cos(varphi_1);
dot_y_B_2 = l_b * dot_varphi_2 * cos(varphi_2);

dot_varphi_B_1 = ((dot_x_B_2 - dot_x_B_1) * cos(varphi_B_2) + (dot_y_B_2 - dot_y_B_1) * sin(varphi_B_2)) / (l_c * sin(varphi_B_1 + varphi_B_2));

dot_x_C = l_b * dot_varphi_1 * sin(varphi_1) + l_c * dot_varphi_B_1 * sin(varphi_B_1);
dot_y_C = l_b * dot_varphi_1 * cos(varphi_1) + l_c * dot_varphi_B_1 * cos(varphi_B_1);

dot_x = [dot_x_C; dot_y_C];
dot_q = [dot_varphi_1; dot_varphi_2];
dot_x = simplify(collect(dot_x, dot_q));

J_qx = simplify(jacobian(dot_x, dot_q));

J = simplify([sin(theta) cos(theta); (cos(theta) / l), -(sin(theta) / l)] * J_qx)
