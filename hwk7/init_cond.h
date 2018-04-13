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


// User Functions
namespace neq {
    const int neqs = 4;
}


namespace init_cond{
    // Define the initial time, timestep, stop time, and the initial RHS array.
    double t, dt, tstop;

    // Define the initial coordinates
    double y_init[neq::neqs]={1.0, 0.0, 0.0, 6.28318530718};

    void prob_init_conds() {
        // Initialize t to zero
        t = 0.0;

        // Initialize stop time to two years
        tstop = 2.0;

        // Ask user for timestep size
        std::cout << "Please input your timestep size (in years): ";
        std::cin >> dt;

        return;
    }
}