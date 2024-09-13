#include <iostream>
#include <thread>
#include <atomic>
#include <string>
#include <webots/Robot.hpp>
#include <webots/Motor.hpp>
#include <webots/PositionSensor.hpp>

#include "leg.hpp"

const robot::FiveLinkParam five_link_param = {
    .l_a = 0.060,
    .l_b = 0.150,
    .l_c = 0.220,
};

const pid::PidConfig pid_leg_length_config = {
    .kp = 2000,
    .ki = 2,
    .kd = 1000000,
    .max_iout = 10,
    .max_out = 100,
};

std::atomic<bool> keep_running(true);
std::atomic<float> l;
std::atomic<float> theta_l;

void read_from_terminal() {
    std::string input;
    while (keep_running) {
        std::cout << "l, theta_l: ";
        std::getline(std::cin, input);

        if (input == "exit") {
            keep_running = false;
            break;
        } else {
            std::istringstream iss(input);
            float temp[2];
            if (iss >> temp[0] >> temp[1]) {
                l = temp[0];
                theta_l = temp[1];
            } else {
                std::cout << "Input Error...\n";
            }
        }
    }
}

int main(int argc, char **argv) {
    auto *robot = new webots::Robot();
    int timeStep = (int)robot->getBasicTimeStep();

    webots::Motor *motor[2] = {
        robot->getMotor("joint_motor_1"),
        robot->getMotor("joint_motor_2"),
    };

    webots::PositionSensor *encoder[2];
    for (int i = 0; i < (int)(sizeof(encoder) / sizeof(encoder[0])); i++) {
        encoder[i] = motor[0]->getPositionSensor();
        encoder[i]->enable(timeStep);
    }
    robot->step(timeStep);
    
    auto *leg = new robot::Leg(five_link_param, 
                               2.8798f - encoder[0]->getValue(),
                               2.8798f - encoder[1]->getValue(),
                               2.8798f - encoder[0]->getValue(),
                               2.8798f - encoder[1]->getValue(),
                               pid_leg_length_config);

    /*std::thread inputThread(read_from_terminal);*/

    while (robot->step(timeStep) != -1 && keep_running) {
        leg->update(2.8798f - encoder[0]->getValue(),
                    2.8798f - encoder[0]->getValue(),
                    0.0f,
                    0.0f);
        leg->ctrl(0.0f, 0.0f, 0.2f);
        motor[0]->setTorque(-leg->left.tau_1);
        motor[1]->setTorque(leg->left.tau_2);
        /*std::cout << "tau_1: " << leg->left.tau_1 << ", tau_2: " << leg->left.tau_2 << ", F_n: " << leg->left.F_n << ", F_l: " << leg->F_ctrl(1) << std::endl;*/
        /*std::cout << "l: " << leg->left.l << std::endl;*/
    };

    keep_running = false;
    /*inputThread.join();*/

    delete leg;
    delete robot;
    return 0;
}
