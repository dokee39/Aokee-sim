#include "leg.hpp"

using namespace std;

namespace robot {
FiveLink::FiveLink(const FiveLinkParam& param,
                   const float varphi_1_init,
                   const float varphi_2_init):
    FiveLinkParam(param),
    varphi_1(varphi_1_init),
    varphi_2(varphi_2_init) {
    forward_solve();
}

void FiveLink::forward_solve() {
    float x_B_1(l_a - l_b * cosf(varphi_1));
    float x_B_2(-l_a + l_b * cosf(varphi_2));
    float y_B_1(l_b * sinf(varphi_1));
    float y_B_2(l_b * sinf(varphi_2));

    float x_B(x_B_1 - x_B_2);
    float y_B(y_B_1 - y_B_2);
    float varphi_B_1(acosf(sqrtf(x_B * x_B + y_B * y_B) / (2 * l_c)) - atan2f(y_B, x_B));

    float x_C(l_a - l_b * cosf(varphi_1) - l_c * cosf(varphi_B_1)); 
    float y_C(l_b * sinf(varphi_1) + l_c * sinf(varphi_B_1));

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
}

void FiveLink::force_virtual_to_actual(const float F_n, const float tau_j) {
}
}
