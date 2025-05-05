function dV = n_last_dV(n,r_ini,r_fin,ToF,v_ini,dv_vector,pars)

%{
###########################################################################
                      Required Last Lambert's Impulse
                        Author: Yago Castilla Lamas
###########################################################################

This function allows for the computation of the last impulse of a
trajectory of n impulses for a given set of the first (n-1) impulses. This
last impulse is computed such that the trajectory satisfies the original
lambert's problem:

INPUTS:
n [scalar] ---> Number of impulses
r_ini [1x3] ---> Initial position of the arc
r_fin [1x3] ---> Final position of the arc
ToF [scalar] ---> Time of flight of the arc
v_ini [1x3] ---> Initial velocity of the arc
dv_vector [(n-1)x3] ---> vector containing the first (n-1) impulses
pars ---> Structure containing problem constants and parameters

OUTPUTS:
dV --> Required last impulse

###########################################################################
Versions Updates:
V2 ---> Modify the impulses to place them at the center of each segment

%}

%Time step
delta_t = ToF/n;

%Current position and velocity
r_vec = r_ini;
v_vec = v_ini;

%Loop on the arcs
for k=1:(n-1)
    %Propagate the orbit for half of a time step
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);

    %Compute velocity after the impulse
    v_vec = v_vec + dv_vector(k,:);

    %Propagate the orbit for half of a time step
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);
end

%Propagate the orbit for half of the last segment (up until last impulse)
[r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);

%Compute the required final impulse
[~,~,~,~,VI] = lambertMR(r_vec, r_fin, delta_t/2, pars.mu_sun, 0, 0, 0, 1);
dV = VI-v_vec;


