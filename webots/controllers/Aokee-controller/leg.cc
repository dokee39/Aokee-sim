#include "leg.hpp"
#include "src/Core/Matrix.h"

using namespace std;

namespace robot {
FiveLink::FiveLink(const FiveLinkParam& param,
                   const float varphi_1_init,
                   const float varphi_2_init):
    FiveLinkParam(param),
    varphi_1(varphi_1_init),
    varphi_2(varphi_2_init) {
    forward_solve();
    jacobian_calc();
}

void FiveLink::forward_solve() {
    x_B_1 = l_a - l_b * cosf(varphi_1);
    x_B_2 = -l_a + l_b * cosf(varphi_2);
    y_B_1 = l_b * sinf(varphi_1);
    y_B_2 = l_b * sinf(varphi_2);

    x_B = x_B_1 - x_B_2;
    y_B = y_B_1 - y_B_2;

    varphi_B_1 = acosf(sqrtf(x_B * x_B + y_B * y_B) / (2 * l_c)) - atan2f(y_B, x_B);

    x_C = l_a - l_b * cosf(varphi_1) - l_c * cosf(varphi_B_1); 
    y_C = l_b * sinf(varphi_1) + l_c * sinf(varphi_B_1);

    l = sqrtf(x_C * x_C + y_C * y_C);
    theta_l = atan2f(x_C, y_C);
}

void FiveLink::inverse_solve() {
    float l_1(sqrtf(l_a * l_a + l * l - 2 * l_a * l * sinf(theta_l))); 
    float l_2(sqrtf(l_a * l_a + l * l + 2 * l_a * l * sinf(theta_l))); 

    varphi_1 = acosf((l_1 * l_1 + l_a * l_a - l * l) / (2 * l_1 * l_a)) + acosf((l_1 * l_1 + l_b * l_b - l_c * l_c) / (2 * l_1 * l_b));
    varphi_2 = acosf((l_2 * l_2 + l_a * l_a - l * l) / (2 * l_2 * l_a)) + acosf((l_2 * l_2 + l_b * l_b - l_c * l_c) / (2 * l_2 * l_b));
}

void FiveLink::jacobian_calc() {
    float l_b_over_l = l_b / l;
    float varphi_B_2 = atan2f(y_C - y_B_2, x_C - x_B_2);
    float sin_varphi_B_1_plus_varphi_B_2 = sinf(varphi_B_1 + varphi_B_2);

    jacobian_trans << -l_b * cosf(theta_l + varphi_B_2) * sinf(varphi_1 - varphi_B_1) / sin_varphi_B_1_plus_varphi_B_2,
                      l_b_over_l * sinf(theta_l + varphi_B_2) * sinf(varphi_1 - varphi_B_1) / sin_varphi_B_1_plus_varphi_B_2,
                      -l_b * cosf(theta_l - varphi_B_1) * sinf(varphi_2 - varphi_B_2) / sin_varphi_B_1_plus_varphi_B_2,
                      l_b_over_l * sinf(theta_l - varphi_B_1) * sinf(varphi_2 - varphi_B_2) / sin_varphi_B_1_plus_varphi_B_2;
}

void FiveLink::force_virtual_to_actual() {
    Eigen::Vector<float, 2> force_virtual;
    Eigen::Vector<float, 2> force_actual;
    force_virtual << F_n, tau_j;
    force_actual = jacobian_trans * force_virtual;
    tau_1 = force_actual(0, 0);
    tau_2 = force_actual(1, 0);
}

void Leg::update(const float left_varphi_1,
                const float left_varphi_2,
                const float right_varphi_1,
                const float right_varphi_2) {
    left.varphi_1 = left_varphi_1;
    left.varphi_2 = left_varphi_2;
    left.forward_solve();
    left.jacobian_calc();
    /*right.varphi_1 = right_varphi_1;*/
    /*right.varphi_2 = right_varphi_2;*/
    /*right.forward_solve();*/
    /*right.jacobian_calc();*/
}

void Leg::ctrl(const float tau_j_l,
               const float tau_j_r,
               const float l_d) {
    F_ctrl(0) = 0.0f;
    F_ctrl(1) = pid_length.calc(left.l, l_d);
    F_ctrl(2) = 100.0f;
    F_ctrl(3) = 0.0f;
    
    F_n = T * F_ctrl;

    left.F_n = F_n(0);
    left.tau_j = tau_j_l;
    left.force_virtual_to_actual();
    right.F_n = F_n(1);
    right.tau_j = tau_j_r;
    right.force_virtual_to_actual();
}
}
