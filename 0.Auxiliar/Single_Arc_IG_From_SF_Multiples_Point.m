function guess = Single_Arc_IG_From_SF_Multiples_Point(r_ini,v_ini,ToF,impulses,m0,nPoints,pars)

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

%Check that the number of points is even
if mod(nPoints,2) ~= 0
    error("Input an even number of points per segment for GPOPS initial guess generation.")
end

%Number of impulses or segments
n_i = length(impulses.times);

%Time of flight of the segments
ToF_segment = ToF/n_i;

%Times of the nodes and initial guess points
guess.times = zeros((nPoints-1)*n_i+1,1);
nodeTimes = linspace(0,ToF,n_i+1);
for k=1:n_i
    for l=1:nPoints-1
    guess.times((nPoints-1)*(k-1)+l) = nodeTimes(k) + (ToF_segment/(nPoints-1))*(l-1);
    end
end
guess.times((nPoints-1)*n_i+1) = nodeTimes(n_i+1);

%Variable to store the state variables at each node (pos,vel)
guess.state = zeros((nPoints-1)*n_i+1,6);

%Variable to store the mass at each node
guess.mass = zeros((nPoints-1)*n_i+1,1);

%Variable to store the thrust at each node
guess.thrust = zeros((nPoints-1)*n_i+1,3);

%Compute initial position, velocity and mass of the arc
r_vec = r_ini;
v_vec = v_ini;
mass = m0;

%Loop on the segments
for k=1:n_i

    %Store the state at the initial node of the segment
    guess.state((nPoints-1)*(k-1)+1,:) = [r_vec,v_vec];

    %Compute the state on the points before the impulse
    for l=1:(nPoints-2)/2
        [r_point,v_point] = propagate_kepler(r_vec,v_vec,ToF_segment*l/(nPoints-1),pars);
        guess.state((nPoints-1)*(k-1)+1+l,:) = [r_point,v_point];
    end

    %Store the mass on this nodes
    for l=1:nPoints/2
        guess.mass((nPoints-1)*(k-1)+l) = mass;
    end

    %Compute the required constant thrust for the segment
    segThrust = impulses.vector(k,:)*1000*mass/ToF_segment;

    %Assign it to all the points of the segment (hold)
    for l=1:nPoints-1
        guess.thrust((nPoints-1)*(k-1)+l,:) = segThrust;
    end

    %Propagate the orbit for the first half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

    %Apply the impulse
    v_vec = v_vec + impulses.vector(k,:);

    %Compute the state on the points after the impulse
    for l=1:(nPoints-2)/2
        [r_point,v_point] = propagate_kepler(r_vec,v_vec,ToF_segment*l/(nPoints-1),pars);
        guess.state((nPoints-1)*(k-1)+nPoints/2+l,:) = [r_point,v_point];
    end
    
    %Propagate the orbit for the second half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

    %Compute the new mass after the impulse
    mass = mass*exp(-norm(impulses.vector(k,:))*1000/(pars.SC.Isp*pars.g0));
    if mass < pars.SC.m0 - pars.SC.mprop
        mass = pars.SC.m0 - pars.SC.mprop;
    end

    %Store the mass of the points after the impulse
    for l=1:(nPoints-2)/2
        guess.mass((nPoints-1)*(k-1)+nPoints/2+l) = mass;
    end

end

%Store last point
guess.state((nPoints-1)*n_i+1,:) = [r_vec,v_vec];
guess.mass((nPoints-1)*n_i+1) = mass;
guess.thrust((nPoints-1)*n_i+1,:) = segThrust;

end
