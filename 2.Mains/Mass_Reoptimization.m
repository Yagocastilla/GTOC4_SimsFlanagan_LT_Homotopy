%{
###########################################################################
                            Mass_Reoptimization
                        Author: Yago Castilla Lamas
               (yago.castilla-lamas@student.isae-supaero.fr)
                              Date: 04/05/2025
###########################################################################

This program loads the results of the optimization performed by the
"Full_Trajectory_Analysis" script. This script reoptimize these results to
maximize the final mass of the S/C. To do this, the script sequentially
reoptimize intervals of adjacent legs of the trajectory. Each interval of
legs is called an optimization window. A constraint is imposed so that the
final velocity of each optimization window is equal to that of the solution
obtained by the "Full_Trajectory_Analysis" script. This is done so that the
solution of the "Full_Trajectory_Analysis" script remains a valid initial
guess for subsequent optimization windows. The script returns a MATLAB
stucture with the reoptimized Sims-Flanagan transcription.

INPUT DATA:
The input data to this script is a MATLAB structure which must be named
"sims_flanagan_trajectory", which contains the Sims-Flanagan transcription
of a given trajectory. The simplest way to use this script is to use
as input data the output file of the "Full_Trajectory_Analysis" script,
which already have the required structure. Moreover, the initial guess used
by the "Full_Trajectory_Analysis" script to generate the trajectory (MATLAB
structure called "Data") must also be provided. To see more details about 
structure of this input data, refer to the description provided in the INPUT
DATA and OUTPUT sections of the "Full_Trajectory_Analysis" script
description.

ALGORITHM CONTROL PARAMETERS:
Several variables are available in the "Algorithm Options" section of the
script to provide the script with the appropiate input data and control its
behaviour:

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
                          "Full_Trajectory_Analysis" script. The precise
                          structure that this file must have is the one
                          described in the OUTPUT part of the "Full_
                          Trajectory_Analysis" script description. Type:
                          string.

resultsAdress ---> Folder adress in which the output file is to be stored.
                   Type: string.

resultsFileName ---> Name of the output file. Type: string.

n_simul ---> Size of the optimization window. Number of arcs that are
             simultaneously optimized. Size: scalar. (recommended value: 2)

OUTPUT:
The output file of this script is the MATLAB structure called "sims_
flanagan_transcription". It is a reoptimized Sims-Flanagan transcription.
Therefore, its structure is exactly the same as the structure of the output
file of the "Full_Trajectory_Analysis" script (see the description of this
script for a more detailed description of this structure).
%}

clear; close all; clc;

%% Algorithm Options

%Select GTOC Solution
solution_index = 4;

%Adress of the file containing the user's initial guess data (only used if
%"solution_index" has been set to 6)
userInitGuess = "../1.GTOC_Data/TrajectoryInitGuess/EDep12_Tour_10.mat";

%File with the optimization results
%(the structure called sims_flanagan_trajectory)
optimization_results = "../3.Results/Original_Moscow/Original_Moscow.mat";

%Adress to store the results and file name
resultsAdress = "../3.Results/Original_Moscow";
resultsFileName = "Original_Moscow_MassReopt.mat";

%Number of arcs simultaneously reoptimized
n_simul = 2;

%% Load GTOC data and trajectory

addpath("../0.Auxiliar")
addpath("../0.Auxiliar/Basic_2_Body")
addpath("../0.Auxiliar/Transcription_Specifics")
addpath("../0.Auxiliar/Transformations")
addpath("../0.Auxiliar/Mass_Optimization")

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
    r_flybys(asteroid,:) = CartesianState(Asteroid_OE, dates(asteroid), pars);
end

%Compute the time of flight of the different arcs
ToFs(1) = (solution.Flyby_Dates(1) - solution.Launch_Date)*86400;
for arc_i=2:length(ToFs)
    ToFs(arc_i) = dates(arc_i)-dates(arc_i-1);
end

%% Mass Reoptimization

%Begin chrono
tic;

%Store the old sims flanagan result
old_sims_flanagan_trajectory = sims_flanagan_trajectory;

%Compute initial position and velocity of the trajectory
r_init = r_earth;
if isfield(sims_flanagan_trajectory, 'vInf')
    v_init = v_earth + sims_flanagan_trajectory.vInf;
else
    v_init = v_launch;
end

%Store the original velocities at each flyby
flyby_vel = zeros(n_arcs,3);
r_vec = r_init;
v_vec = v_init;
for arc_i = 1:n_arcs
    %Propagate the arc to obtain final position and velocity of the arc
    [r_vec, v_vec] = propagate_segmented(r_vec,v_vec,sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).vector,sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).times,ToFs(arc_i),pars);
    
    %Store the flyby velocity
    flyby_vel(arc_i,:) = v_vec;
end

%Variable to store the impulses
impulses = zeros(n_i,3,n_arcs);

%Transform the initial guess into a more suitable form
for arc_i=1:n_arcs
    impulses(:,:,arc_i) = sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).vector;
end

%Variable to store the masses
arc_mass = zeros(1,n_arcs+1);

%Initial conditions
r_vec = r_init;
v_vec = v_init;
arc_mass(1) = pars.SC.m0;

%Optimization loop with simultaneous optimization
for arc_i=1:n_arcs-n_simul+1

    %Compute the parameters of the optimization window
    r_ini = r_vec; %Initial position of the optimization window
    v_ini = v_vec; %Initial velocity of the optimization window
    r_obj = r_flybys(arc_i:arc_i+n_simul-1,:); %Extract the position of the flybys in the optimization window
    arc_ToFs = ToFs(arc_i:arc_i+n_simul-1); %Extract the time of flight of the arcs in the optimization window
    initMass = arc_mass(arc_i); %Extract the SC masses at each arc in the optimization window
    comb_guess = impulses(:,:,arc_i:arc_i+n_simul-1); %Initial guess for the impulses

    %Reoptimization of the optimization window
    combined_sims_flanagan = Mass_Simultaneous_NLP(r_ini,v_ini,r_obj,flyby_vel(arc_i+n_simul-1,:),arc_ToFs,initMass,comb_guess,n_i,pars.SC.T,pars);

    %Store the results
    impulses(:,:,arc_i:arc_i+n_simul-1) = combined_sims_flanagan.vector;
    times(:,arc_i:arc_i+n_simul-1) = combined_sims_flanagan.times;

    %Store the unchanged results
    sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).vector = impulses(:,:,arc_i);
    sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).times = times(:,arc_i);

    %Propagate the arc
    [r_vec, v_vec] = propagate_segmented(r_vec,v_vec,impulses(:,:,arc_i),times(:,arc_i),ToFs(arc_i),pars);

    %Update the masses
    for k=1:n_simul
        arc_mass(arc_i+k) = updateArcMass(combined_sims_flanagan.vector(:,:,k).*1000,arc_mass(arc_i),pars);
        if arc_mass(arc_i+k) < pars.SC.m0 - pars.SC.mprop
            arc_mass(arc_i+k) = pars.SC.m0 - pars.SC.mprop;
        end
    end

    %Display progress message
    fprintf("Progress. Arc %d finished. \n", arc_i);
end

%Stop chrono and store the reoptimization time
executionTime = toc;
sims_flanagan_trajectory.reoptTime = executionTime;

%% Save the results
SaveAdress = fullfile(resultsAdress, resultsFileName);
save(SaveAdress, "sims_flanagan_trajectory");
