#include <iostream>
#include <thread>
#include <atomic>
#include <string>
#include <webots/Robot.hpp>
#include <webots/Motor.hpp>
#include <webots/PositionSensor.hpp>

#include "leg.hpp"

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
    
    const robot::FiveLinkParam five_link_param = {
        .l_a = 60,
        .l_b = 150,
        .l_c = 220,
    };
    auto *five_link = new robot::FiveLink(five_link_param, 2.8798f - encoder[0]->getValue(), 2.8798f - encoder[1]->getValue());
    l = five_link->l;
    theta_l = five_link->theta_l;

    std::thread inputThread(read_from_terminal);

    while (robot->step(timeStep) != -1 && keep_running) {
        /*std::cout << five_link->l << ", " << five_link->theta_l << std::endl;*/
        five_link->l = l;
        five_link->theta_l = theta_l;
        five_link->inverse_solve();

        motor[0]->setPosition(2.8798f - five_link->varphi_1);
        motor[1]->setPosition(five_link->varphi_2 - 2.8798f);

        for (int i = 0; i < 200 / timeStep; i++) {
            robot->step(timeStep);
        }
    };

    keep_running = false;
    inputThread.join();

    delete five_link;
    delete robot;
    return 0;
}
