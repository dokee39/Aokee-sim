#include <iostream>
#include <webots/Robot.hpp>
#include <webots/Motor.hpp>
#include <webots/PositionSensor.hpp>

#include "leg.hpp"

int main(int argc, char **argv) {
    webots::Robot *robot = new webots::Robot();
    webots::Motor *motor[2] = {
        robot->getMotor("joint_motor_1"),
        robot->getMotor("joint_motor_2"),
    };
    webots::PositionSensor *encoder[2];
    for (int i = 0; i < (int)(sizeof(encoder) / sizeof(encoder[0])); i++) {
        encoder[i] = motor[0]->getPositionSensor();
        encoder[i]->enable(1);
    }

    /*const robot::FiveLinkParam five_link_param = {*/
    /*    .l_a = 60,*/
    /*    .l_b = 150,*/
    /*    .l_c = 220,*/
    /*};*/
    /*robot::FiveLink five_link(five_link_param, 2.8798f - encoder[0]->getValue(), 2.8798f + encoder[1]->getValue());*/

    int timeStep = (int)robot->getBasicTimeStep();
    while (robot->step(timeStep) != -1) {
        std::cout << "varphi_1 = " << 2.8798f - encoder[0]->getValue() << ", "\
                  << "varphi_2 = " << 2.8798f - encoder[0]->getValue() << std::endl;
        motor[0]->setTorque(0);
        motor[1]->setTorque(0);
        for (int i = 0; i < 200 / timeStep; i++) {
            robot->step(timeStep);
        }
    };

    delete robot;
    return 0;
}
