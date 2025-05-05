function dV_lamb = next_lambert(r_ini,v_ini,dv_ini,r_obj,ToF,ToF_segment,pars)

%{
###########################################################################
                      Required Next Lambert's Impulse
                        Author: Yago Castilla Lamas
###########################################################################

This function apply to a two impulses manoeuvre that satisfies a lambert's
problem. The function takes the initial impulse as an input and returns the
required value of the second impulse.

INPUTS:
r_ini ---> Initial position vector [km]
v_ini ---> Initial velocity vector [km/s]
dv_ini ---> Initial impulse (vector) [km/s]
r_obj ---> Final position vector of the lambert's arc [km]
ToF ---> Time of flight of the lambert arc [s]
ToF_segment ---> Time of flight of the first segment [s]
pars ---> Structure containing problem data including mu

OUTPUT:
dV_lamb ---> Magnitude of the second impulse [km/s]

NOTE: to avoid errors while executing "lambertMR" function, this function
computes the eccentricty of the orbit after the first impulse. If the orbit
is hyperbolic, the program exits the function and the manoeuvre is assigned
a prohibitive value. 
%}

%Provide the impulse
r_vec = r_ini;
v_vec = v_ini + dv_ini;

%Check if the orbit is hyperbolic
[~, e, ~, ~, ~, ~] = Build_OE(r_vec, v_vec, pars.mu_sun);
if e>1
    dV_lamb = 10e4;
    return
end

%Propagate the orbit up to the next impulse
[r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment,pars);

%Compute the remaining time
Remaining_ToF = ToF - ToF_segment;

%Compute the Lambert Arc to reach the final position
[~,~,~,~,VI] = lambertMR(r_vec, r_obj, Remaining_ToF, pars.mu_sun, 0, 0, 0, 1);

%Compute the next lambert's impulse
dV_lamb = norm(VI - v_vec);

