function [r_fin,v_fin] = propagate_kepler(r_ini,v_ini,delta_t,pars)

%{
###########################################################################
                      Required Next Lambert's Impulse
                        Author: Yago Castilla Lamas
###########################################################################

This program propagates the orbit of an object for a specified amount
of time:

INPUTS:
r_ini ---> Initial position vector [km]
v_ini ---> Initial velocity vector [km/s]
delta_t ---> Time interval of propagation [s]
pars ---> Structure containing problem data including mu

OUTPUT:
r_fin ---> Final position vector [km]
v_fin ---> Final velocity vector [km/s]
%}

%Obtain orbital elements
[a, e, i, RAN, Omega, MA_i] = Build_OE(r_ini, v_ini, pars.mu_sun);
Body_OE = [a/pars.AU, e, i, RAN, Omega, MA_i];

%Propagate orbit
[r_fin, v_fin] = CartesianState(Body_OE, delta_t, pars);

