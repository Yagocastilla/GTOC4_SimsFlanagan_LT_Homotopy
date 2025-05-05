function guess = Single_Arc_IG_From_SF(r_ini,v_ini,ToF,impulses,m0,pars)

%{
###########################################################################
              Single Arc Generate GPOPS Initial Guess From SF
                        Author:  Yago Castilla Lamas
###########################################################################

This program takes a sims-flanagan transcription of an arc, represented by
the impulses magnitude and direction and the times where these take place,
and converts it into a suitable structure to be used by GPOPS.

In this specific method, each impulse is converted to an equivalent thrust
vector that is considered to be applied constantly along each segment. The
value of the thrust is provided in a zero order hold form. That is, up
until the next node, the thrust is considered constant and equal to the
value provided in the initial node of the current segment.

The trajectory, on the other hand, is computed applying the instataneous
impulses and propagating the trajectory as a keplerian motion. The mass is
updated after each segment using tschiolkovsky's equation and this is the
one used to compute the thrust of the next segment.

INPUTS:
r_ini ---> Initial position vector [km]
v_ini ---> Initial velocity vector [km/s]
ToF ---> Time of flight of the arc [s].
impulses ---> Structure containing two informations about the impulses
              (output of the single arcs SF setups functions):
              impulses.vector: Matrix comtaining the impulses of the arc.  
                               The rows are the impulses and the columns   
                               the cartesian components. Size Nx3. [km/s]
              impulses.times: Times relative to the beginning of the       
                              segment where the impulses take place. Size 
                              Nx1. [s]
m0 ---> Initial mass of the arc [Kg]
pars ---> Structure containing problem constants and parameters.

OUTPUT:
guess ---> Structure containing the initial guess discretized in the nodes
           between the segments of the SF transcription provided through
           the "impulses" structure (guess.times, guess.state, guess.mass,
           guess.thrust).

###########################################################################
Versions Updates:
V2 ---> -Modify the impulses to place them at the center of each segment
        -Zero order hold for the thrust
        -Can only be used with a uniform discretization in time for the
        nodes in the arc.

%}

%Number of impulses or segments
n_i = length(impulses.times);

%Time of flight of the segments
ToF_segment = ToF/n_i;

%Structure for the initial guess
guess.times = zeros(2*n_i,1);
nodeTimes = linspace(0,ToF,n_i+1);
for k=1:n_i
    guess.times(2*k-1) = nodeTimes(k);
    guess.times(2*k) = nodeTimes(k+1) - 1;
end
guess.times(2*n_i) = nodeTimes(n_i+1);

%Variable to store the state variables at each node (pos,vel)
guess.state = zeros(2*n_i,6);

%Variable to store the mass at each node
guess.mass = zeros(2*n_i,1);

%Variable to store the thrust at each node
guess.thrust = zeros(2*n_i,3);

%Compute initial position, velocity and mass of the arc
r_vec = r_ini;
v_vec = v_ini;
mass = m0;

%Loop on the segments
for k=1:n_i

    %Store the state at the initial node of the segment
    guess.state(2*k-1,:) = [r_vec,v_vec];

    %Store the mass at the initial node
    guess.mass(2*k-1) = mass;

    %Compute the required constant thrust for the segment
    segThrust = impulses.vector(k,:)*1000*mass/ToF_segment;

    %Assign it to both the initial and final nodes of the segment (hold)
    guess.thrust(2*k-1,:) = segThrust;
    guess.thrust(2*k,:) = segThrust;

    %Propagate the orbit for the first half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

    %Apply the impulse
    v_vec = v_vec + impulses.vector(k,:);
    
    %Propagate the orbit for the second half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

    %Store the position at the final node of the segment
    guess.state(2*k,:) = [r_vec,v_vec];

    %Compute the new mass after the segment
    mass = mass*exp(-norm(impulses.vector(k,:))*1000/(pars.SC.Isp*pars.g0));
    if mass < pars.SC.m0 - pars.SC.mprop
        mass = pars.SC.m0 - pars.SC.mprop;
    end

    %Store the mass at the final node of the segment
    guess.mass(2*k) = mass;

end

end
