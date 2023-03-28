This is an application to generate geometry and calculate volume for Miura and Kresling Bellows 
at different stable states. This application consists of 4 steps；

Step 1: Choose n and plot design space
Users need to input an integer n (n>=3). The geometric design space is shown on the upper right 
figure according to Reid et al.
(https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.013002). 

Step 2: Select a crease pattern
Kresling pattern is the default pattern, d1/R=2pi/n, users do not need to input d1/R. If switch 
to Miura pattern, users need to input d1/R by using the slider or do it manually.

Step 3: Compute number of stable states
Users need to enter the rest design parameters (\phi1,\phi2,R,m), and then the numer of stable 
states will be indicated in the box. If the structure is flat-foldable, it will be indicated.

Step 4: Output geometry and compute volume
This step is to plot the strcutures at all feasible stable states, and will also output the height 
ratio, radial expansion ratio, volume expansion ratio etc.


