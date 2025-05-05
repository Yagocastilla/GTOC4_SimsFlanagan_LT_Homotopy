function [r_vec, v_vec] = propagate_segmented(r_ini,v_ini,dv_vector,times,ToF,pars)

%{
###########################################################################
                          Propagate Segmented Arc
                        Author:  Yago Castilla Lamas
###########################################################################

This function takes an arc discretized in N segments and provides the final
position and velocity of the S/C afer the arc.

INPUTS:
r_ini ---> Initial position vector [km]
v_ini ---> Initial velocity vector [km/s]
impulses ---> Matrix comtaining the impulses of the arc. The rows are the
              impulses and the columns the cartesian components. Size Nx3.
times ---> Times relative to the beginning of the segment where the
           impulses take place. Size Nx1.
ToF ---> Time of flight of the arc [s].

OUTPUTS:
r_vec ---> Final position vector [km]
v_vec ---> Final velocity vector [km/s]

###########################################################################
Versions Updates:
V2 ---> It can handle arcs where the first impulse is not placed at t=0

%}

%Number of segments
n_segments = size(dv_vector,1);

%Current position and velocity
r_vec = r_ini;
v_vec = v_ini;

%Include the time of flight in the matrix of times
times = [times;ToF];

%Propagate the orbit until the first impulse
[r_vec,v_vec] = propagate_kepler(r_vec,v_vec,times(1),pars);

%Loop on the segments
for k=1:n_segments
    %Compute velocity after the impulse
    v_vec = v_vec + dv_vector(k,:);

    %Compute the time of flight of the segment
    ToF_segment = times(k+1)-times(k);

    %Propagate the orbit for a time step
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment,pars);
end

