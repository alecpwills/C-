/*
File: ode_rhs.h
Author: Alec Wills
Date: April 11 2018
Editor: VSCode
Purpose: To allow access to the ode_rhs() function, that the user can update according to
         their system's needs.
*/


// Includes
#include<cmath>



// User Functions
namespace ode_rhs{
    double* ode_rhs(double t, double dt, double* y0){        
        // Define return array
        double* rhs = new double[neq::neqs];
        
        // Define Gravitational constant and M_sun in solar mass/AU (Keplerian) units
        const double GRAV = 39.47, MSUN = 1.0;

        // Define and initialize the initial condition parameters to update,
        // unpacking the y0 input array
        // Here
        double x = y0[0];
        double y = y0[1];
        double vx = y0[2];
        double vy = y0[3];

        // Intialize the updated variables
        rhs[0] = vx; // RHS of x-coordinate equation, the incremental change in x
        rhs[1] = vy; // RHS of y-coordinate equation, the incremental change in y
        rhs[2] = -GRAV*MSUN*x/std::pow( (std::pow(x,2) + std::pow(y, 2)), 1.5); // RHS of x velocity equation, increment
        rhs[3] = -GRAV*MSUN*y/std::pow( (std::pow(x,2) + std::pow(y, 2)), 1.5); // RHS of y velocity equation, increment
        
        // Return updated RHS
        return rhs;
    }
}

