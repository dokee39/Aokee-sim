#pragma once

#include <Eigen>
#include <cmath>

namespace robot {
struct FiveLinkParam {
    float l_a;
    float l_b;
    float l_c;
};

class FiveLink: public FiveLinkParam {
public:
    explicit FiveLink(const FiveLinkParam& param,
                      const float varphi_1_init,
                      const float varphi_2_init);
    ~FiveLink() = default;

    void forward_solve();
    void inverse_solve();
    void jacobian_calc();
    void force_virtual_to_actual(const float F_n, const float tau_j);

    float l;
    float theta_l;
    float varphi_1;
    float varphi_2;
    float tau_1;
    float tau_2;

private:
    Eigen::Matrix<float, 2, 2> jacobian;
};

class Leg {
public:
    explicit Leg(const FiveLinkParam& five_link_param,
                 const float leg_left_varphi_1,
                 const float leg_left_varphi_2,
                 const float leg_right_varphi_1,
                 const float leg_right_varphi_2):
        leg_left(five_link_param, leg_left_varphi_1, leg_left_varphi_2),
        leg_right(five_link_param, leg_right_varphi_1, leg_right_varphi_2) {
    }
    ~Leg();

    FiveLink leg_left;
    FiveLink leg_right;
};
}

