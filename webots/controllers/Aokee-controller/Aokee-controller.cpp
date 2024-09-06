#include <webots/Robot.hpp>
#include <webots/Motor.hpp>

int main(int argc, char **argv) {
    webots::Robot *robot = new webots::Robot();
    webots::Motor *motor[2] = {
        robot->getMotor("joint_motor_1"),
        robot->getMotor("joint_motor_2"),
    };

    int timeStep = (int)robot->getBasicTimeStep();
    while (robot->step(timeStep) != -1) {
        motor[0]->setPosition(0.5);
        motor[1]->setPosition(-0.5);
    };

    delete robot;
    return 0;
}
