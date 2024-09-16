syms theta_w_1(t) theta_w_2(t) theta_l_1(t) theta_l_2(t) theta_b(t)
syms dot_theta_w_1 dot_theta_w_2 dot_theta_l_1 dot_theta_l_2 dot_theta_b
syms ddot_theta_w_1 ddot_theta_w_2 ddot_theta_l_1 ddot_theta_l_2 ddot_theta_b

syms R_w R_l l_1 l_2 l_w_1 l_w_2 l_b_1 l_b_2 l_c g m_w m_l m_b I_w I_l_1 I_l_2 I_b I_z

diffs = [diff(theta_w_1, t, t) diff(theta_w_2, t, t) diff(theta_l_1, t, t) diff(theta_l_2, t, t) diff(theta_b, t, t) ...
         diff(theta_w_1, t) diff(theta_w_2, t) diff(theta_l_1, t) diff(theta_l_2, t) diff(theta_b, t)];
dots = [ddot_theta_w_1 ddot_theta_w_2 ddot_theta_l_1 ddot_theta_l_2 ddot_theta_b ...
        dot_theta_w_1 dot_theta_w_2 dot_theta_l_1 dot_theta_l_2 dot_theta_b];
square_of_dots = [dot_theta_w_1^2 dot_theta_w_2^2 dot_theta_l_1^2 dot_theta_l_2^2 dot_theta_b^2];
zeros = [0 0 0 0 0];

s_w = R_w / 2 *(theta_w_1 + theta_w_2);
ddot_s_w = simplify(subs(subs(diff(s_w, t, t), diffs, dots), square_of_dots, zeros));

h_b = l_1 / 2 * cos(theta_l_1) + l_2 / 2 * cos(theta_l_2);
ddot_h_b = simplify(subs(subs(diff(h_b, t, t), diffs, dots), square_of_dots, zeros));

s_b = R_w / 2 * (theta_w_1 + theta_w_2) + l_1 / 2 * sin(theta_l_1) + l_2 / 2 * sin(theta_l_2);
ddot_s_b = simplify(subs(subs(diff(s_b, t, t), diffs, dots), square_of_dots, zeros));

s_l_1 = R_w * theta_w_1 + l_w_1 * sin(theta_l_1);
ddot_s_l_1 = simplify(subs(subs(diff(s_l_1, t, t), diffs, dots), square_of_dots, zeros));

s_l_2 = R_w * theta_w_2 + l_w_2 * sin(theta_l_2);
ddot_s_l_2 = simplify(subs(subs(diff(s_l_2, t, t), diffs, dots), square_of_dots, zeros));

h_l_1 = h_b - l_b_1 * cos(theta_l_1);
ddot_h_l_1 = simplify(subs(subs(diff(h_l_1, t, t), diffs, dots), square_of_dots, zeros));

h_l_2 = h_b - l_b_2 * cos(theta_l_2);
ddot_h_l_2 = simplify(subs(subs(diff(h_l_2, t, t), diffs, dots), square_of_dots, zeros));

phi = R_w / (2 * R_l) * (- theta_w_1 + theta_w_2) - l_1 / (2 * R_l) * sin(theta_l_1) + l_2 / (2 * R_l) * sin(theta_l_2);
ddot_phi = simplify(subs(subs(diff(phi, t, t), diffs, dots), square_of_dots, zeros));


syms tau_w_1(t) tau_w_2(t) tau_j_1(t) tau_j_2(t)
syms F_ws_1 F_ws_2 F_wh_1 F_wh_2 F_ns_1 F_ns_2 F_nh_1 F_nh_2 f_1 f_2

vars = [ddot_theta_w_1 ddot_theta_w_2 ddot_theta_l_1 ddot_theta_l_2 ddot_theta_b ...
        tau_w_1 tau_w_2 tau_j_1 tau_j_2];
binding_force = [F_ws_1 F_ws_2 F_wh_1 F_wh_2 F_ns_1 F_ns_2 F_nh_1 F_nh_2 f_1 f_2];
ddot_thetas = [ddot_theta_w_1 ddot_theta_w_2 ddot_theta_l_1 ddot_theta_l_2 ddot_theta_b];

eqn1 = m_w * R_w * ddot_theta_w_1 == f_1 - F_ws_1;
eqn2 = m_w * R_w * ddot_theta_w_2 == f_2 - F_ws_2;

eqn3 = I_w * ddot_theta_w_1 == tau_w_1 - f_1;
eqn4 = I_w * ddot_theta_w_2 == tau_w_2 - f_2;

eqn5 = m_l * ddot_s_l_1 == F_ws_1 - F_ns_1;
eqn6 = m_l * ddot_s_l_2 == F_ws_2 - F_ns_2;

eqn7 = m_l * ddot_h_l_1 == F_wh_1 - F_nh_1 - m_l * g;
eqn8 = m_l * ddot_h_l_2 == F_wh_2 - F_nh_2 - m_l * g;

eqn9 = I_l_1 * ddot_theta_l_1 == - tau_w_1 + tau_j_1 - (F_ws_1 * l_w_1 + F_ns_1 * l_b_1) * cos(theta_l_1) + (F_wh_1 * l_w_1 + F_nh_1 * l_b_1) * sin(theta_l_1);
eqn10 = I_l_2 * ddot_theta_l_2 == - tau_w_2 + tau_j_2 - (F_ws_2 * l_w_2 + F_ns_2 * l_b_2) * cos(theta_l_2) + (F_wh_2 * l_w_2 + F_nh_2 * l_b_2) * sin(theta_l_2);

eqn11 = m_b * ddot_s_b == F_ns_1 + F_ns_2;

eqn12 = m_b * ddot_h_b == F_nh_1 + F_nh_2 - m_b * g;

eqn13 = I_b * ddot_theta_b == - tau_j_1 - tau_j_2 - (F_ns_1 + F_ns_2) * l_c * cos(theta_b) + (F_nh_1 + F_nh_2) * l_c * sin(theta_b);

eqn14 = I_z * ddot_phi == (- f_1 + f_2) * R_l;

eqn15 = F_wh_1 == F_wh_2;

thetas = [sin(theta_l_1) sin(theta_l_2) sin(theta_b) ...
          cos(theta_l_1) cos(theta_l_2) cos(theta_b)];
% thetas_alt = [theta_l_1(t) theta_l_2(t) theta_b(t) (1 - theta_l_1^2) (1 - theta_l_2^2) (1 - theta_b^2)];
thetas_alt = [theta_l_1(t) theta_l_2(t) theta_b(t) 1 1 1];
% thetas_alt = [0 0 0 1 1 1];
thetas_mul = [theta_l_1^2 theta_l_2^2 theta_l_1*theta_l_2 ...
              theta_w_1^2 theta_l_2^2 theta_w_1*theta_w_2];
zeros = [0 0 0 0 0 0];

eqns_binding_force = [eqn1 eqn2 eqn3 eqn4 eqn5 eqn6 eqn7 eqn8 eqn12 eqn15];
eqns_binding_force = subs(eqns_binding_force, thetas, thetas_alt);
eqns = [eqn9 eqn10 eqn11 eqn13 eqn14];
eqns = subs(eqns, thetas, thetas_alt);

sol = solve(eqns_binding_force, binding_force);
disp(sol);
sols = simplify(collect([sol.F_ws_1 sol.F_ws_2 sol.F_wh_1 sol.F_wh_2 sol.F_ns_1 sol.F_ns_2 sol.F_nh_1 sol.F_nh_2 sol.f_1 sol.f_2], vars));

results = simplify(collect(subs(eqns, binding_force, sols), vars));
% disp(results);
results = simplify(subs(results, thetas_mul, zeros));
% disp(results);

results = solve(results, ddot_thetas);
disp(results);

