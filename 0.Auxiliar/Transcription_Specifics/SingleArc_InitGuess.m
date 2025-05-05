function arc_guess = SingleArc_InitGuess(r_ini,v_ini,r_fin,ToF,n_i,pars)

%{
###########################################################################
                         Single Arc Initial Guess
                        Author: Yago Castilla Lamas
###########################################################################

This programm generates an initial guess for the optimization of an arc
using a Sims Flanagan transcription of n_i segments. Four different methods
to generate the initial guess can be selected. All the methods are based in
solving one or several lambert's problems. More details on each method can
be found on the report of this work. The initial guess does not necessarily
satisfy all the constraints of the arc (final position, maximum impulse,
etc).

INPUTS:
r_ini ---> Initial position vector [km]
v_ini ---> Initial velocity vector [km/s]
r_fin ---> Final position vector [km]
ToF ---> Time of flight of the arc [s]
n_i ---> Number of impulses or segments
pars ---> Structure containing problem constants and parameters
init ---> Flag defining the method to generate the initial guess:
          1 - Divided lambert's impulse
          2 - Succesive lambert's impulses
          3 - Scaled succesive lambert's impulses
          4 - Minimize next lambert's impulse
maxdV ---> Maximum impulse that the propulsion system can provide in a
           segment (only needed when using initial guess 4) [km/s]

OUTPUTS:
arc_guess.vector ---> Matrix (n_i x 3). The rows are the different impulses
                      and the columns correspond to the 3 cartesian
                      components.
arc_guess.times ---> Vector containing the times at wich the impulses are
                     performed. These times are measured taking as t=0 the
                     beginning of the arc.

%}

%Variable to store the initial condition
dv_0 = zeros(n_i,3);

%Lambert's arc required initial velocity
[~,~,~,~,VI_lamb] = lambertMR(r_ini, r_fin, ToF, pars.mu_sun, 0, 0, 0, 1);

%Required lambert's impulse
lamb_impulse = VI_lamb - v_ini;

r_vec = r_ini; %Initial position vector
v_vec = v_ini; %Initial velocity vector

ToF_segment = ToF/n_i; %Time of flight of the segment

%Loop on the segments
for k=1:size(dv_0,1)
    %Propagate the trajectory for half of the segment
    %(up to where the impulse is placed)
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

    %Remaining time of flight
    Remaining_ToF = ToF - (2*k-1)*ToF_segment/2;

    %Lambert's arc required initial velocity
    [~,~,~,~,VI_lamb_segment] = lambertMR(r_vec, r_fin, Remaining_ToF, ...
                                          pars.mu_sun, 0, 0, 0, 1);

    %Required lambert's impulse
    lamb_impulse_segment = VI_lamb_segment - v_vec;

    %Scale the impulse and set it as initial condition
    dv_0(k,:) = lamb_impulse_segment/(n_i-k+1);

    %Provide the impulse
    v_vec = v_vec + dv_0(k,:);

    %Propagate the trajectory for the remainig half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);
end

%Store the results
arc_guess.vector = dv_0;

%Times of the impulses
arc_guess.times = zeros(size(arc_guess.vector,1),1);
for j=1:length(arc_guess.times)
    arc_guess.times(j) = (2*j-1)*ToF/(2*n_i);
end

end
