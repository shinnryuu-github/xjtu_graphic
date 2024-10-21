#include "solver.h"

#include <Eigen/Core>

using Eigen::Vector3f;

// External Force does not changed.

// Function to calculate the derivative of KineticState
KineticState derivative(const KineticState& state)
{
    return KineticState(state.velocity, state.acceleration, Eigen::Vector3f(0, 0, 0));
}

// Function to perform a single Forward Euler step
KineticState forward_euler_step([[maybe_unused]] const KineticState& previous,
                                const KineticState& current)
{
    extern float time_step;
    
    Eigen::Vector3f new_position = current.position + current.velocity * time_step;
    Eigen::Vector3f new_velocity = current.velocity + current.acceleration * time_step;

    return KineticState(new_position, new_velocity, current.acceleration);
}

// Function to perform a single Runge-Kutta step
KineticState runge_kutta_step([[maybe_unused]] const KineticState& previous,
                              const KineticState& current)
{
    extern float time_step;

    KineticState k1 = current;
    KineticState k2 = KineticState(current.position + k1.velocity * time_step / 2,
                                              current.velocity + k1.acceleration * time_step / 2,
                                              current.acceleration);
    KineticState k3 = KineticState(current.position + k2.velocity * time_step / 2,
                                              current.velocity + k2.acceleration * time_step / 2,
                                              current.acceleration);
    KineticState k4 = KineticState(current.position + k3.velocity * time_step,
                                              current.velocity + k3.acceleration * time_step,
                                              current.acceleration);

    Eigen::Vector3f new_position = current.position + (k1.velocity + 2 * k2.velocity + 2 * k3.velocity + k4.velocity) * (time_step / 6);
    Eigen::Vector3f new_velocity = current.velocity + (k1.acceleration + 2 * k2.acceleration +2 * k3.acceleration + k4.acceleration) * (time_step / 6);

    return KineticState(new_position, new_velocity, current.acceleration);
}

// Function to perform a single Backward Euler step
KineticState backward_euler_step([[maybe_unused]] const KineticState& previous,
                                 const KineticState& current)
{
    extern float time_step;

    Eigen::Vector3f new_velocity = current.velocity + current.acceleration * time_step;
    Eigen::Vector3f new_position = current.position + new_velocity * time_step;

    return KineticState(new_position, new_velocity, current.acceleration);
}

// Function to perform a single Symplectic Euler step
KineticState symplectic_euler_step(const KineticState& previous, const KineticState& current)
{
    (void)previous;
    extern float time_step;

    Eigen::Vector3f new_velocity = current.velocity + current.acceleration * time_step;
    Eigen::Vector3f new_position = current.position + new_velocity * time_step;

    return KineticState(new_position, new_velocity, current.acceleration);
}
