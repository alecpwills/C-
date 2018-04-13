/*
File: ode.cpp
Author: Alec Wills
Date: April 11 2018
Editor: VSCode
Purpose: To integrate ODE coupled equations using either RK2 or RK4 methods, and write out
         the results to a file for plotting.
*/


// Includes
#include<iostream>
#include<fstream>
#include "init_cond.h"
#include "ode_rhs.h"
/*
int main(){
    init_cond::prob_init_conds();
    double* arr = init_cond::y_init;
    std::cout << init_cond::t << init_cond::tstop << init_cond::dt << arr[0] << arr[1] << arr[2]<<arr[3]<<std::endl;

    std::cout << init_cond::y_init[0] << init_cond::y_init[1] << init_cond::y_init[2] << init_cond::y_init[3];
    return 0;
} */

// Function Skeletons

double* add_arr(double*, double*);
double* scalar_mult(double, double*);
void rk2();
void rk4();

int main() {
    // Initialize the conditions specified in header
    std::cout << "Setting initial conditions...\n";
    init_cond::prob_init_conds();
    std::cout << "Initial conditions set.\n";
    std::cout << "You may view these conditions in the init_cond.h header file and update as needed for the given problem.\n";
    std::cout << "\n" << "The user has specified " << neq::neqs << " coupled equations.\n";
    std::cout << "You may view the incrementing functions in ode_rhs.h header, within the C++ function there.\n";
    std::cout << "These functions are incremented by dt when stepping forward in the Euler method,\n and are here in array form to represent the function-vector F(y^n, t) in the Runge-Kutta algorithm.\n";
    std::cout << "\n";
    // Ask user to choose between RK2 and RK4
    int choice;
    std::cout << "Please choose an integration algorithm:\n";
    std::cout << "1: RK2 (Second Order Runge-Kutta)\n";
    std::cout << "2: RK4 (Fourth Order Runge-Kutta)\n";
    std::cin >> choice;
    std::cout << "\n";

    switch (choice){
        case 1: rk2();
        case 2: rk4();
    }
    

return 0;
}


// User Functions

double* add_arr(double* arr1, double* arr2){
    int len = sizeof(arr1)/sizeof(arr1[0]);
    double* returnarr = new double[len];
    for (int i=0; i<len; i++){
        returnarr[i] = arr1[i]+arr2[i];
    }
    return returnarr;
}

double* scalar_mult(double scalar, double* arr){
    // Length of array
    int len = sizeof(arr)/sizeof(arr[0]);

    // Return array
    double* newarr = new double[len];

    // Loop over indices and multiply
    for (int i=0; i < len; i++){
        newarr[i] = scalar*arr[i];
    }
    return newarr;
}

void rk2(){
    // Will start with the initial condition array and
    // step through the evolution, writing the new RHS
    // at each step

    // Define the y array that will be updated throughout
    double* y = init_cond::y_init;
    
    // Define the t values to input into ode_rhs
    double t = init_cond::t, tstop = init_cond::tstop;
    double dt = init_cond::dt;

    // Open file to write to
    std::ofstream rk2file;
    rk2file.open("rk2_output.dat");
    std::cout << "Writing position to rk2_output.dat: t x y format." << std::endl;
    while (t < tstop){
        // Generate the RHS for input y
        // that are incremented by dt
        double* RHS = ode_rhs::ode_rhs(t, dt, y);

        // Define the yn that will be combined with k1
        double y1[neq::neqs];

        // Scalar multiply the RHS by the given dt for k1
        double k1[neq::neqs];
        for (int i=0; i<neq::neqs; i++){
            k1[i] = dt*RHS[i];
        }

        // Delete old RHS so that it can be reused in next loop, and since the function returns
        // a dynamically assigned array
        delete [] RHS;

        // Initialize the y1 vector as defined
        for (int i=0; i<neq::neqs; i++){
            y1[i] = y[i] + 0.5*k1[i];
        }
        // Reevaulate the RHS with the updated y1 vector
        double* RHSnew = ode_rhs::ode_rhs(t+0.5*dt, dt, y1);

        // Update y array values
        // Total reassignment not allowed, loop over elements
        for (int i=0; i<neq::neqs; i++){
            y[i] = y[i] + dt*RHSnew[i];
        }

        // Delete old RHS new so that it can be used in next loop
        delete [] RHSnew;

        // Update t value
        t = t + dt;

        // Write out to the console, can be redirected to file
        rk2file << t << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << "\n";
    }
    std::cout << "Closing file.\n";
    rk2file.close();
    std::cout << "File closed. Integration complete."<< std::endl;
    return; 
} 

void rk4(){
    // Will start with the initial condition array and
    // step through the evolution, writing the new RHS
    // at each step

    // Define the y array that will be updated throughout
    double* y = init_cond::y_init;
    
    // Define the t values to input into ode_rhs
    double t = init_cond::t, tstop = init_cond::tstop;
    double dt = init_cond::dt;

    // Open file to write to
    std::ofstream rk4file;
    rk4file.open("rk4_output.dat");
    std::cout << "Writing position to rk4_output.dat: t x y format." << std::endl;
    while (t < tstop){
        // Generate the RHS for input y
        double* RHS = ode_rhs::ode_rhs(t, dt, y);

        // Define the necessary yn that will be used as
        // input arguments for updating the RHS
        double y2[neq::neqs], y3[neq::neqs], y4[neq::neqs];

        // Scalar multiply RHS array by the given dt for k1
        double k1[neq::neqs];
        for (int i=0; i<neq::neqs; i++){
            k1[i] = dt*RHS[i];
        }

        // Delete old RHS so that it can be reused in next loop
        delete [] RHS;

        // Initialize y2 to be used in the RHS argument
        // for finding k2
        for (int i=0; i<neq::neqs; i++){
            y2[i] = y[i] + 0.5*k1[i];
        }

        // Reevaulate the RHS for use in finding k2
        double* RHS2 = ode_rhs::ode_rhs(t+0.5*dt, dt, y2);

        // Step 2, evaluate k2 as dt*RHS2
        double k2[neq::neqs];
        for (int i=0; i<neq::neqs; i++){
            k2[i] = dt*RHS2[i];
        }

        // Delete RHS2 to be updated next loop
        delete [] RHS2;

        // Initialize y3 to be used as argument in k3
        for (int i=0; i<neq::neqs; i++){
            y3[i] = y[i] + 0.5*k2[i];
        }

        // Reevaluate the RHS for use in finding k3
        double* RHS3 = ode_rhs::ode_rhs(t+0.5*dt, dt, y3);

        // Evaluate k3 at dt*RHS3
        double k3[neq::neqs];
        for (int i=0; i<neq::neqs; i++){
            k3[i] = dt*RHS3[i];
        }

        // Delete RHS3 to be updated next loop
        delete [] RHS3;

        // Initialize y4 to be used as argument in k4
        for (int i=0; i<neq::neqs; i++){
            y4[i] = y[i] + k3[i];
        }

        // Reevaluate the RHS for use in finding k4
        double* RHS4 = ode_rhs::ode_rhs(t+dt, dt, y4);

        // Evaluate k4 as dt*RHS4
        double k4[neq::neqs];
        for (int i=0; i<neq::neqs; i++){
            k4[i] = dt*RHS4[i];
        }

        // Update y array values
        // Total reassignment not allowed, loop over elements
        for (int i=0; i<neq::neqs; i++){
            y[i] = y[i] + (1.0/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
        }

        // Delete old RHS new so that it can be used in next loop
        delete [] RHS4;

        // Update t value
        t = t + dt;

        // Write out to the console, can be redirected to file
        rk4file << t << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << "\n";
    }
    std::cout << "Closing file.\n";
    rk4file.close();
    std::cout << "File closed. Integration complete."<< std::endl;
    return; 
} 