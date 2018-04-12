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
#include "neq.h"
#include "init_cond.h"



// User Functions

namespace ode_rhs{
    double ode_rhs(double t, double dt, double y0[neq::neq]){
        // Define return array
        double rhs[neq::neq];
        
        // Define Gravitational constant and M_sun in solar mass/AU (Keplerian) units
        const double GRAV = 39.47d0, MSUN = 1.0d0;

        // Define and initialize the initial condition parameters to update,
        // unpacking the y0 input array
        double x = y0[0];
        double y = y0[1];
        double vx = y0[2];
        double vy = y0[3];

        // Intialize the updated variables
        rhs[0] = x + dt*vx; // RHS of x-coordinate equation
        rhs[1] = y + dt*vy; // RHS of y-coordinate equation
        rhs[2] = vx -dt*GRAV*MSUN*x/std::pow(x**2+y**2, 1.5d0) // RHS of x velocity equation
        rhs[3] = vy -dt*GRAV*MSUN*y/std::pow(x**2+y**2, 1.5d0) // RHS of y velocity equation
        
        // Return updated RHS
        return rhs;
    }
}

