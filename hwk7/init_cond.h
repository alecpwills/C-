/*
File: init_cond.h
Author: Alec Wills
Date: April 11 2018
Editor: VSCode
Purpose: To allow the user to input specified intial conditions, regardless of the form of
         the coupled ODE system. It allows for a variable number of coupled equations, which
         determines the size of the initial condition array.
*/


// Includes
#include<iostream>
#include "neq.h"



// User Functions
namespace init_cond{
    // Define the initial time, timestep, stop time, and the initial RHS array.
    double t, dt, tstop, y_init[neq::neq];

    void prob_init_conds() {
        // Initialize t to zero
        t = 0.0d0

        // Initialize stop time to two years
        tstop = 2.0d0

        // Ask user for timestep size
        std::cout << "Please input your timestep size (in years): ";
        std::cin >> dt;

        // Initialize the initial x-coordinate, in AU
        y_init[0] = 1.0d0

        // Initialize the initial y-coordinate, in AU
        y_init[1] = 0.0d0

        // Initialize the initial x velocity
        y_init[2] = 0.0d0

        // Initialize the initial y velocity
        y_init[3] = 6.28318530718d0  // AU/year, i.e. the circumference of the Earth's orbit/year

        return;
    }
}