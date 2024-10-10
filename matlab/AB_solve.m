function [A, B] = AB_solve(A_sym, B_sym, param_num)
    syms R_w R_l l_1 l_2 l_w_1 l_w_2 l_b_1 l_b_2 l_c g m_w m_l m_b I_w I_l_1 I_l_2 I_b I_z
    
    param = [R_w R_l l_1 l_2 l_w_1 l_w_2 l_b_1 l_b_2 l_c g m_w m_l m_b I_w I_l_1 I_l_2 I_b I_z];

    A = double(subs(A_sym, param, param_num));
    B = double(subs(B_sym, param, param_num));
end
