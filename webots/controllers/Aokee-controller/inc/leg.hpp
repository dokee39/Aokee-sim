#pragma once

#include <Eigen>
#include <cmath>
#include "pid_controller.hpp"
#include "src/Core/Matrix.h"

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
    void force_virtual_to_actual();

    float l;
    float theta_l;
    float varphi_1;
    float varphi_2;
    float tau_1 = 0.0f;
    float tau_2 = 0.0f;
    float F_n = 0.0f;
    float tau_j = 0.0f;

private:
    float varphi_B_1;
    float x_B_1;
    float x_B_2;
    float y_B_1;
    float y_B_2;
    float x_B;
    float y_B;
    float x_C;
    float y_C;
    Eigen::Matrix<float, 2, 2> jacobian_trans;
};

class Leg {
public:
    explicit Leg(const FiveLinkParam& five_link_param,
                 const float left_varphi_1,
                 const float left_varphi_2,
                 const float right_varphi_1,
                 const float right_varphi_2,
                 const pid::PidConfig& pid_length_config):
        left(five_link_param, left_varphi_1, left_varphi_2),
        right(five_link_param, right_varphi_1, right_varphi_2),
        pid_length(pid_length_config) {
        T << 1, 1, 1, -1, -1, 1, 1, 1;
    }
    ~Leg() = default;

    void update(const float left_varphi_1,
                const float left_varphi_2,
                const float right_varphi_1,
                const float right_varphi_2);
    void ctrl(const float tau_j_l,
              const float tau_j_r,
              const float l_d);

    FiveLink left;
    FiveLink right;

    pid::Pid pid_length;

    // [F_n,l F_n,r]^T
    Eigen::Vector<float, 2> F_n;
    // [F_psi F_l F_gravity F_inertial]^T
    Eigen::Vector<float, 4> F_ctrl;
    
    Eigen::Matrix<float, 2, 4> T;
};
}

