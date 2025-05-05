%{
###########################################################################
                      ResultsAnalysis_GPOPSInitGuess
                        Author: Yago Castilla Lamas
               (yago.castilla-lamas@student.isae-supaero.fr)
                              Date: 04/05/2025
###########################################################################

This program loads the results of the "Full_Trajectory_Analysis" or
"Mass_Reoptimization" scripts to analyze the results and provide the user
with figures that allows him to quickly asses them. Moreover, it generates
a structure which contains an initial guess to be used on a higher fidelity
optimization tool. This initial guess is a sequence of discretized points
of the trajectory. On each point, the values of the state vector, the mass
and the thrust are provided.

The outputs of the script are then: a structure with the S/C control
history; a graph with the mean thrust of each segment of the trajectory; a
graph showing the evolution of the SC mass; a 3D view of the trajectory;
and a MATLAB structure which contains the initial guess for the higher
fidelity solver. The script also checks that the provided Sims-Flanagan
transcription perform all the flybys of the different asteroids in time and
therefore is a valid trajectory. Finally, the script computes the
difference in velocity between the SC and the final asteroid, which should
be as small as possible so that the rendez-vous condition is satisfied.

INPUT DATA:
The input data to this script is a MATLAB structure which must be named
"sims_flanagan_trajectory" which contains the Sims-Flanagan transcription
of a given trajectory. The simplest way to use this script is to use
as input data the output file of the "Full_Trajectory_Analysis" or
"Mass_Reoptimization" scripts, which already have the required structure.
Moreover, the initial guess used by the "Full_Trajectory_Analysis" script
to generate the trajectory (MATLAB structure called "Data") must also be
provided. To see more details about the structure of this input data refer
to the description provided in the INPUT DATA and OUTPUT sections of the
"Full_Trajectory_Analysis" script description.

This data is provided through three variables defined in the "Input Files"
section of this script. These three variables are:

solution_index ---> Flag to select the solution under study. Size: scalar.
                    Five example solutions are already implemented:
                    1 - Moscow
                    2 - Barbee
                    3 - Johnson
                    4 - Original Moscow
                    5 - Impulsive Johson
                    If the user is studying a self-defined trajectory,
                    this parameter should be set to 6 and the adress of the
                    structure containing the initial guess must be defined
                    in the "userInitGuess" parameter.

userInitGuess ---> Adress of the file containing input MATLAB structure
                   containing the data of initial guess defined by the
                   user, used by the "Full_Trajectory_Analysis" script to
                   generate the optimized trajectory. This parameter is
                   used only if the solution_index parameter has been set
                   to 6. The precise structure that this file must have is
                   the one described in the INPUT DATA part of the
                   "Full_Trajectory_Analysis" script description. Type:
                   string.

optimization_results ---> Adress of the file with the MATLAB structure
                          containing the results of the optimization of the
                          Sims-Flanagan transcription performed by the
                          "Full_Trajectory_Analysis" or "Mass_Reoptimization"
                          scripts. The precise structure that this file
                          must have is the one described in the OUTPUT part
                          of the "Full_Trajectory_Analysis" script
                          description. Type: string.

OUTPUT:
There is several ouput structures for this script which allows the user to
evaluate the optimized trajectory:

control_history ---> A structure containing all the information about the
                     impulses to be performed by the spacecraft and some
                     additional information about the Sims-Flanagan segment
                     corresponding to that impulse:
                        
                     control_history.impulse ---> The cartesian components
                     of the impulses to be applied to the S/C. Size:
                     n_segments x 3 (n_segments = n_i x n_arcs) [km/s].

                     control_history.impulse_magnitude ---> The magnitude
                     of each of the provided impulses. Size: n_segments x 1
                     [km/s].

                     control_history.date ---> The date at which the
                     impulse is provided. Size: n_segments x 1 [MJD].
            
                     control_history.times ---> The elapsed time since the
                     launch date to the impulse, in seconds. Size:
                     n_segments x 1 [s].

                     control_history.ToF_segments ---> The ToF of the
                     corresponding segment, in seconds. Size: n_segments
                     x 1 [s].

                     control_history.mass ---> The mass of the spacecraft
                     for the corresponding segment. Size: n_segments x 1
                     [kg].

                     control_history.thrust ---> The mean thrust magnitude
                     of the segment computed as the magnitude of the
                     impulse divided by the ToF of the segment and
                     multiplied by the mass of the S/C. Size: n_segment
                     x 1 [N].

sol_guess ---> A structure suitable to be used as initial guess for a
               higher fidelity solver such as GPOPS. This structure
               contains the information about a sequence of discretized
               points all over the trajectory. In particular, each leg of
               the trajectory is discretized in 51 points. Fields:

               sol_guess.state ---> The cartesian components of the
               position and velocity of the SC at each point. The structure
               of each entry is: [x,y,z,dx/dt,dy/dt,dz/dt] in km and km/s.
               Size: 51 x 6 x n_arcs (last index corresponds to the leg of
               the trajectory in which the point is located).

               sol_guess.times ---> The elapsed time since the launch date.
               Size: 51 x 1 x n_arcs [s].

               sol_guess.mass ---> The mass of the S/C at each point. Size:
               51 x 1 x n_arcs [kg].

               sol_guess.thrust ---> The cartesian components of the thrust
               vector which is being applied at each point. Size: 51 x 3
               x n_arcs [N].

               sol_guess.dateModifiers ---> The number of days each flyby
               has been shifted due to the trajectory optimization,
               compared to the initial guess. This result is directly
               copied from the input data. To obtain the optimized flyby
               dates, it suffices with adding this vector to the provided
               initial guess for the flyby dates (solution.Flyby_Dates +
               sol_guess.dateModifiers). Size: 1 x n_arcs [days].

               sol_guess.vInf ---> Cartesian components of the earth
               departure excess velocity. This result is directly copied
               from the input data. To obtain the intial velocity of
               the S/C wrt the Sun, this velocity should be added to the
               Earth velocity at the launch date. Size: 1 x 3 [km/s].

               sol_guess.removedAsteroids ---> The index of the asteroids
               that have been removed from the provided initial guess
               asteroid sequence. This result is directly copied from the
               input data. Size: 1D array.

err_flyby ---> The distance between the S/C and the asteroids on each of
               the flybys. Size: 1 x n_arcs [km].

vRendezvous ---> The magnitude of the difference between the S/C velocity
                 on the last flyby and the velocity of the last asteroid on
                 the date of this flyby. Size: 3 x 1 [km/s].


OUTPUT GRAPHS:
- Mean Thrust Plot:
This graph shows the S/C mean thrust of each segment along the trajectory.
The red dashed lines on this figure represent each of the asteroids flybys.
The blue dashed line represents the maximum thrust that the S/C propulsion
system can provide. For the trajectory to be feasible, the mean thrust of
every segment should be under this line.

- Mass Plot:
This graph shows the S/C mass evolution along the trajectory. It should be
remembered that in the applied transcription, the mass of the S/C is
constant along each leg and it is only updated between legs. The red dashed
lines equally represent the asteroids flybys. The blue dashed line in this
case represents the empty mass of the spacecraft. Therefore, for the
trajectory to be feasible in terms of mass, the mass of the S/C should
remain above this line for the full trajectory.

- Trajectory 3D Plot:
In this graph, a 3D view of the trajectory is shown. Each of the asteroid
flybys is marked with a star.
%}

clear; close all; clc;

%% Input Files

%Select solution to be converted (set to 6 for a new trajectory defined by
%the user)
solution_index = 4;

%Adress of the file containing the user's initial guess data (only used if
%"solution_index" has been set to 6)
userInitGuess = "../1.GTOC_Data/TrajectoryInitGuess/EDep12_Tour_10.mat";

%File with the optimization results
%(the structure called sims_flanagan_trajectory)
optimization_results = "../3.Results/Original_Moscow/Original_Moscow_MassReopt.mat";

%% Load GTOC data and trajectory

addpath("../0.Auxiliar")
addpath("../0.Auxiliar/Basic_2_Body")
addpath("../0.Auxiliar/Transcription_Specifics")
addpath("../0.Auxiliar/Transformations")

%Load GTOC problem data
pars = load("../1.GTOC_Data/pars.mat").pars;

%Load trajectory
if solution_index == 1
    solution_name = "Moscow";
    GTOC_sol = "../1.GTOC_Data/TrajectoryInitGuess/Moscow_Sol.mat";
elseif solution_index == 2
    solution_name = "Barbee";
    GTOC_sol = "../1.GTOC_Data/TrajectoryInitGuess/Barbee_Sol.mat";
elseif solution_index == 3
    solution_name = "Johnson";
    GTOC_sol = "../1.GTOC_Data/TrajectoryInitGuess/Johnson_Sol.mat";
elseif solution_index == 4
    solution_name = "Original Moscow";
    GTOC_sol = "../1.GTOC_Data/TrajectoryInitGuess/Moscow_Original_Sol.mat";
elseif solution_index == 5
    solution_name = "Impulsive_Johnson";
    GTOC_sol = "../1.GTOC_Data/TrajectoryInitGuess/Impulsive_Johnson_Sol.mat";
elseif solution_index == 6
    solution_name = "User";
    GTOC_sol = userInitGuess;
else
    error("Invalid entry for solution_index parameter.")
end
solution = load(GTOC_sol).Data;

%Loading txt file with Asteroid Data
GTOC4Asteroids = readtable("../1.GTOC_Data/GTOC4.txt");
Asteroids_OE = GTOC4Asteroids(solution.Asteroid_Seq_ID,3:8);

%% Load optimization results and check number of arcs

%Load optimization results structure
sims_flanagan_trajectory = load(optimization_results).sims_flanagan_trajectory;

%Remove asteroids
if isfield(sims_flanagan_trajectory, 'removedAsteroids')
    solution = removeAsteroids(solution, sims_flanagan_trajectory.removedAsteroids);
end

%Number of arcs in the solution
n_arcs = length(solution.Flyby_Dates);

%Number of arcs in the optimization results
n_arcs_opt = sum(startsWith(fieldnames(sims_flanagan_trajectory), "Arc"));

%Check that both numbers are equal
if n_arcs ~= n_arcs_opt
    error("The number of arcs in the trajectory and in the " + ...
          "optimization results does not correspond.")
end

%Number of segments per arc
n_i = length(sims_flanagan_trajectory.Arc1.times);

%% Earth Launch 

%Earth orbital elements on reference date (54000 MJD)
OE_Earth_0 = [0.999988049532578, 0.01671681163160, 0.0008854353079654,...
    175.40647696473, 287.61577546182, 257.60683707535];

%Earth position on the launch date
[r_earth,v_earth] = CartesianState(OE_Earth_0, (solution.Launch_Date-pars.Epoch.Earth_Ephem)*86400, pars);

%First asteroid orbital elements on reference date (54800 MJD)
Asteroid1_OE = table2array(Asteroids_OE(1,:));

%Position of the first asteroid on S/C flyby date
r_asteroid = CartesianState(Asteroid1_OE, (solution.Flyby_Dates(1)-pars.Epoch.Asteroid_Ephem)*86400, pars);

%Variable to store the ToF of the arcs
ToFs = zeros(length(solution.Flyby_Dates),1);

%Time of flight of the arc earth-first asteroid
ToFs(1) = (solution.Flyby_Dates(1) - solution.Launch_Date)*86400;

%Compute the inital lambert's arc departure velocity (to be used only if
%there is no optimization of vInf)
[~,~,~,~,v_launch] = lambertMR(r_earth, r_asteroid, ToFs(1), pars.mu_sun, 0, 0, 0, 1);

%% Modify the dates of the solution

if isfield(sims_flanagan_trajectory, 'dateModifiers')
    solution.Flyby_Dates = solution.Flyby_Dates + sims_flanagan_trajectory.dateModifiers;
end

%% Asteroid flybys position vectors and TOFs

%Vector to store the position vectors
r_flybys = zeros(length(solution.Flyby_Dates),3);
v_asteroids = zeros(length(solution.Flyby_Dates),3);

%Vector to store the elapsed time since the reference date for each flyby
dates = zeros(length(solution.Flyby_Dates),1);

%Loop on the asteroids
n_asteroids = length(solution.Flyby_Dates);
for asteroid=1:n_asteroids
    %Asteroid orbital elements on reference date (54800 MJD)
    Asteroid_OE = table2array(Asteroids_OE(asteroid,:));

    %Elapsed time since the reference date until the flyby date in seconds
    dates(asteroid) = (solution.Flyby_Dates(asteroid)-pars.Epoch.Asteroid_Ephem)*86400;
    
    %Position of the asteroid on S/C flyby
    [r_flybys(asteroid,:), v_asteroids(asteroid,:)] = CartesianState(Asteroid_OE, dates(asteroid), pars);
end

%Compute the time of flight of the different arcs
ToFs(1) = (solution.Flyby_Dates(1) - solution.Launch_Date)*86400;
for arc_i=2:length(ToFs)
    ToFs(arc_i) = dates(arc_i)-dates(arc_i-1);
end

%% Obtain control history from the SF transcription

%Number of arcs in the SF results
n_arcs_SF = sum(startsWith(fieldnames(sims_flanagan_trajectory), "Arc"));

%Total number of segments
n_segments = n_arcs_SF*n_i;

%Control history: impulse vectors
control_history.impulse = zeros(n_segments,3);
for arc_i=1:n_arcs
    control_history.impulse((arc_i-1)*n_i+1:arc_i*n_i,:) = sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).vector;
end

%Control history: dates and ToF of the segments
control_history.date = zeros(n_segments,1);

%First arc segments
ToF_segment = ToFs(1)/n_i/86400;
for j=1:n_i
    control_history.ToF_segments(j) = ToF_segment*86400; %[s]
end
for j=1:n_i
    control_history.date(j) = solution.Launch_Date + (2*j-1)*ToF_segment/2;
end

%Rest of the arcs
for arc_i=2:n_arcs_SF
    %Get the flyby dates
    Flyby1_Date = solution.Flyby_Dates(1,arc_i-1); %[MJD]
    Flyby2_Date = solution.Flyby_Dates(1,arc_i); %[MJD]

    %Time of flight of the arc and of the segments
    ToF_arc = Flyby2_Date - Flyby1_Date;
    ToF_segment = ToF_arc/n_i; %[MJD]
    for j=1:n_i
        control_history.ToF_segments((arc_i-1)*n_i+j) = ToF_segment*86400; %[s]
    end

    %Dates of the impulses in MJD
    for j=1:n_i
        control_history.date((arc_i-1)*n_i+j) = Flyby1_Date + (2*j-1)*ToF_segment/2;
    end
end

%Control history: times since the launch date in MJD
control_history.times = control_history.date-solution.Launch_Date*ones(length(control_history.date),1);

%Control history: impulse magnitude
control_history.impulse_magnitude = zeros(n_segments,1);
for j=1:n_segments
    control_history.impulse_magnitude(j) = norm(control_history.impulse(j,:));
end

%Control history: constant mass at each arc
control_history.mass = zeros(n_segments,1);
SCMass = pars.SC.m0;
control_history.mass(1:n_i,1) = SCMass;
for arc_i=2:n_arcs_SF
    SCMass = updateArcMass(control_history.impulse((arc_i-2)*n_i+1:(arc_i-1)*n_i,:).*1000,SCMass,pars);
    if SCMass < pars.SC.m0 - pars.SC.mprop
        SCMass = pars.SC.m0 - pars.SC.mprop;
    end
    control_history.mass((arc_i-1)*n_i+1:arc_i*n_i,1) = SCMass;
end

%Control history: segment mean thrust
control_history.thrust = zeros(n_segments,1);
for j=1:n_segments
    control_history.thrust(j) = control_history.mass(j)*control_history.impulse_magnitude(j)*1000/control_history.ToF_segments(j); %[N]
end

%% Commanded Mean Thrust Plot

ThrustPlotLimit=0.3; % [N]

%Mean thrust over time
figure1 = figure;
axes1 = axes('Parent',figure1); % Create axes
hold(axes1,'on');

plot(control_history.times, control_history.thrust,'LineWidth',2,...
     'Color',[1 0 0],HandleVisibility = "Off") %Plot commanded thrust

%Plot maximum thrust
line([control_history.times(1) control_history.times(end)],...
     [pars.SC.T pars.SC.T], 'Parent',axes1,'LineStyle','--',...
     Color = "b", DisplayName='Maximum Available Thrust')

%Plot asteroid flybys
line([solution.Flyby_Dates(1)-solution.Launch_Date solution.Flyby_Dates(1)-solution.Launch_Date],...
     [0 ThrustPlotLimit],'Parent',axes1,'LineStyle','--', Color = "r",...
     DisplayName='Asteroid Flybys')

for iL=2:length(solution.Flyby_Dates)
line([solution.Flyby_Dates(iL)-solution.Launch_Date solution.Flyby_Dates(iL)-solution.Launch_Date],...
     [0 ThrustPlotLimit],'Parent',axes1,'LineStyle','--', Color = "r",...
     HandleVisibility = "Off")
end

%Axis and legend
ylabel('Thrust [N]');
xlabel('Time [days]');
title(sprintf("%s's Trajectory Commanded Mean Thrust",solution_name))
legend()
ylim(axes1,[0 ThrustPlotLimit]);
xlim(axes1,[0 3600])
box(axes1,'on');
hold(axes1,'off');

%% SC Mass Plot

%Mass over time
figure2 = figure;
axes2 = axes('Parent',figure2); % Create axes
hold(axes2,'on');

plot(control_history.times, control_history.mass,'LineWidth',2,...
     'Color',[1 0 0],DisplayName = "SC Mass") %Plot SC mass

%Plot empty mass
line([control_history.times(1) control_history.times(end)],...
     [pars.SC.m0-pars.SC.mprop pars.SC.m0-pars.SC.mprop], 'Parent',axes2,'LineStyle','--',...
     Color = "b", DisplayName='Empty Mass')

%Plot asteroid flybys
line([solution.Flyby_Dates(1)-solution.Launch_Date solution.Flyby_Dates(1)-solution.Launch_Date],...
     [0 pars.SC.m0],'Parent',axes2,'LineStyle','--', Color = "r",...
     DisplayName='Asteroid Flybys')

for iL=2:length(solution.Flyby_Dates)
line([solution.Flyby_Dates(iL)-solution.Launch_Date solution.Flyby_Dates(iL)-solution.Launch_Date],...
     [0 pars.SC.m0],'Parent',axes2,'LineStyle','--', Color = "r",...
     HandleVisibility = "Off")
end

%Axis and legend
ylabel('Mass [Kg]');
xlabel('Time [days]');
title(sprintf("%s's Trajectory SC Mass",solution_name))
legend()
ylim(axes2,[0 pars.SC.m0]);
xlim(axes2,[0 3600]);
box(axes2,'on');
hold(axes2,'off');

%% Trajectory Computation

%Number of points per segment
mesh_refiner = 100;

%Variable to store the trajectory
segments_xyz = zeros(n_segments,mesh_refiner,3); %[arc,point,coordinate]

%Compute initial position and velocity of the trajectory
r_vec = r_earth;
if isfield(sims_flanagan_trajectory, 'vInf')
    v_vec = v_earth + sims_flanagan_trajectory.vInf;
else
    v_vec = v_launch;
end

%Number of points in the first half of each segment
num_points_first = floor(mesh_refiner/2);

%Variable to store state vector at the end of each leg
legsFinalState = zeros(n_arcs_SF,6);

%Loop on the segments
for k=1:n_segments

    %Time of flight of the segment
    ToF_segment = control_history.ToF_segments(k);

    %Obtain the trajectory of the first half of the segment
    segments_xyz(k,1:num_points_first,:) = obtain_Segment_Trajectory(r_vec,v_vec,ToF_segment/2,...
                                            pars,num_points_first);

    %Propagate the orbit for the first half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

    %Apply the impulse
    v_vec = v_vec + control_history.impulse(k,:);

    %Obtain the trajectory of the second half of the segment
    segments_xyz(k,num_points_first+1:mesh_refiner,:) = obtain_Segment_Trajectory(r_vec,v_vec,ToF_segment/2,...
                                            pars,mesh_refiner-num_points_first);

    %Propagate the orbit for the second half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

    %Store final position and velocity of each leg
    if mod(k, n_i) == 0
        legsFinalState(k/n_i,:) = [r_vec, v_vec];
    end
end

%Store the final velocity of the spacecraft
vF_SC = v_vec;

%Reshape the vector to have a unique column with all the points
trajectory = zeros((n_segments+1)*mesh_refiner,3);
for k=1:n_segments
    trajectory(mesh_refiner*(k-1)+1:mesh_refiner*k,:) = segments_xyz(k,:,:);
end

%% 3D Graphical Representation

figure (3)
hold on
plot3(trajectory(:,1).',trajectory(:,2).',trajectory(:,3).', 'k-', DisplayName='Trajectory'); % Trajectory
plot3(r_earth(1), r_earth(2), r_earth(3), 'r.', 'MarkerSize', 8, DisplayName='Earth'); % Initial Position
plot3(r_flybys(end,1), r_flybys(end,2), r_flybys(end,3), 'rpentagram', 'MarkerSize', 8, DisplayName='Final Asteroid'); % Final Position

plot3(r_flybys(1,1), r_flybys(1,2), r_flybys(1,3), 'k*', 'MarkerSize', 8, DisplayName='Asteroids'); % First Asteroid
for asteroid=2:n_asteroids-1
    plot3(r_flybys(asteroid,1), r_flybys(asteroid,2), r_flybys(asteroid,3),...
          'k*', 'MarkerSize', 8, HandleVisibility = "Off"); % Asteroids
end

plot3(0, 0, 0, 'bo', 'MarkerSize', 10, DisplayName='Sun'); % Sun Representation

%Axes and title
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
xlim([-2.5e8 2.5e8]);
axis equal
grid on;
title(sprintf("%s's Trajectory",solution_name));
legend('FontSize', 6)
hold off

%% Generate initial guess structure for GPOPS

%Number of points per segment of GPOPS initial guess
nPointsGPOPS = 6;

%Create the structure of the initial guess for the complete trajectory
sol_guess.times = zeros((nPointsGPOPS-1)*n_i+1,1,n_arcs);
sol_guess.state = zeros((nPointsGPOPS-1)*n_i+1,6,n_arcs);
sol_guess.mass = zeros((nPointsGPOPS-1)*n_i+1,1,n_arcs);
sol_guess.thrust = zeros((nPointsGPOPS-1)*n_i+1,3,n_arcs);

%Initial position, velocity and mass
r_vec = r_earth;
if isfield(sims_flanagan_trajectory, 'vInf')
    v_vec = v_earth + sims_flanagan_trajectory.vInf;
else
    v_vec = v_launch;
end
SCMass = pars.SC.m0;

%Loop on the arcs
for arc_i=1:n_arcs

    %Time of flight of the arc
    ToF = ToFs(arc_i);

    %Get the initial mass of the arc
    m0 = SCMass;

    %Create the structure of the initial guess
    guess = Single_Arc_IG_From_SF_Multiples_Point(r_vec,v_vec,ToF,sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)),m0,nPointsGPOPS,pars);

    %Insert the obtained guess in the complete solution structure
    sol_guess.times(:,1,arc_i) = guess.times;
    sol_guess.state(:,:,arc_i) = guess.state;
    sol_guess.mass(:,1,arc_i) = guess.mass;
    sol_guess.thrust(:,:,arc_i) = guess.thrust;

    %Propagate the arc to obtain initial position and velocity of next arc
    [r_vec, v_vec] = propagate_segmented(r_vec,v_vec,sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).vector,sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).times,ToF,pars);

    %Update the mass
    SCMass = updateArcMass(sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).vector.*1000,SCMass,pars);
    if SCMass < pars.SC.m0 - pars.SC.mprop
        SCMass = pars.SC.m0 - pars.SC.mprop;
    end
end

%If the dates have been modified, store it also in the guess structure
if isfield(sims_flanagan_trajectory, 'dateModifiers')
    sol_guess.dateModifiers = sims_flanagan_trajectory.dateModifiers;
end

%If the v infinity have been optimized, store it also in the guess
if isfield(sims_flanagan_trajectory, 'vInf')
    sol_guess.vInf = sims_flanagan_trajectory.vInf;
end

%If asteroids have been removed, store it also in the guess structure
if isfield(sims_flanagan_trajectory, 'removedAsteroids')
    sol_guess.removedAsteroids = sims_flanagan_trajectory.removedAsteroids;
end

%% Check that the trajectory satisfy the flybys

%Compute the error on the position of each flyby
err_flyby = zeros(1,n_arcs);
for arc_i=1:n_arcs
    %Final position of the arc
    seg = n_i*arc_i;
    final_r = [segments_xyz(seg,end,1), segments_xyz(seg,end,2), segments_xyz(seg,end,3)];

    %Flyby error
    err_flyby(arc_i) = norm(final_r - r_flybys(arc_i,:));
end

%Compute velocity of last asteroid
Asteroid_OE = table2array(Asteroids_OE(end,:));
[r_endAsteroid, v_endAsteroid] = CartesianState(Asteroid_OE, (solution.Flyby_Dates(end)-pars.Epoch.Asteroid_Ephem)*86400, pars);

%Compute the final velocity difference between S/C and final asteroid
vRendezvous = norm(v_endAsteroid - vF_SC);
fprintf("The difference in velocity at final rendezvous is: %d \n",vRendezvous)

%Check the that the error in the flybys is under 10 km
error_flag = 0;
for arc_i=1:n_arcs
    if err_flyby(arc_i)>10
        error_flag = 1;
    end
end
if error_flag==1
    error("The trajectory does not satisfy the flybys. Arc %d \n",arc_i)   
else
    disp("The trajectory satisfy the flybys.")
end
