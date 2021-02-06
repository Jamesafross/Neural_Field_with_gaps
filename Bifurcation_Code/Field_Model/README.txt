To compute the bifurcation curves for the spatially extended model, XPPAUTO was used to solve and continue the solutions to the characteristic equation and implicit function theorem. To reproduce these curves, you will need to run each of the files provided in XPP. XPP will not simulate the equations as they are not actually differential equations. 

Use the Monte Carlo method to find the fixed point (solution of the set of equations), by selecting 'Sing pts'->'monte(C)ar'. XPP will ask for parameters for the Monte Carlo search and search windows for the model parameters. The solutions are output into the XPP terminal. 

Next change the initial conditions to the found solution and press go. This will 'simulate' the equations for these initial condition.

To find the solution as a function of one of the parameters, open AUTO ('File'->'Auto'). AUTO opens in a new window. You will need change the 'Numerics' to suitable values. Then press 'Run' and the solution will be continued as a function of the first parameter. You can change which parameter you continue in by changing the 'Axes'.