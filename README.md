# Gillespie-SIS
SIS Epidemic model simulation with Gillespie's algorithm (MATLAB)

These files were created for academic research performed under the supervision of Dr. Eric Forgoston at Montclair State University, Department of Mathematical Sciences.

The goal of said research was to explore the theory of "Early Warning Signals" and its application to mathematical epidemiology and population models with Allee effect, thereby allowing for the study of control mechanisms to affect the behavior of these dynamical systems. 

The model chosen and implemented (in both MATLAB and FORTRAN) in the code found in this repository follow the SIS (susceptible infected susceptible) epidemic model.

The simulations found herein were accomplished by employing Gillespie's algorithm to simulate stochastic epidemiological data.

-----------------------------------------------------------------------------------------------------------------------------
Gillespie_SIS_V3 is the main simulation file. It plots the results of the simulation and allows the user to "pulse" the population at a given point.

Gillespie_SIS_V6 is similar to V3, but only used to generate a predetermined number of time series.

Gillespie_SIS_V7 is similar to V3, however, the pulse occurs when autocorrelation values reach a statistical threshold. This is accomplished by reading from .csv files containing said values (imported from R).

Gillespie_SIS_V5 is a FORTRAN version of V3 without the implementaion of control measures.
