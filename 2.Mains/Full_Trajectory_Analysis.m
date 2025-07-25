%{
###########################################################################
                          Full Trajectory Analysis
                        Author: Yago Castilla Lamas
               (yago.castilla-lamas@student.isae-supaero.fr)
                              Date: 04/05/2025
###########################################################################

This script is designed to perform the optimization of asteroid tours
trajectories designed to solve the problem proposed by GTOC4. The script
takes as input the ordered list of the n_asteroids asteroids to be visited,
as well as the launch date and an initial guess for the flyby dates of each
of the asteroids, which should be sufficiently good. The script then
performs the optimization of the low-thrust trajectory by applying the
Sims-Flanagan transcription. The complete methodology is described in
detail in the report that comes with this tool. However, a brief
description of the workflow is provided here. The script works by
sequentially optimizing intervals of adjacent legs of the trajectory. Each
interval of legs is called an optimization window. For each optimization
window, the following steps are performed:

1 - The flyby dates of the remaining legs of the trajectory are reoptimized
    using the 1-ISF transcription.

2 - The 1-ISF transcription is used to obtain an optimal final velocity for
    the optimization window. This velocity is obtained so that the
    feasibility of subsequent legs is maximized.

3 - The multi-impulse Sims-Flanagan transcription of the optimization
    window is optimized to minimize the difference between the actual final
    velocity of the optimization window and the optimal computed one while
    respecting all the problem constraints such as the maximum thrust.

The script returns a MATLAB stucture with the optimized Sims-Flanagan
trasncription of the trajectory which also contains the information
required to obtained the optimized flyby dates obtained by the algorithm.

INPUT DATA:
The input data to this script is an initial guess for the asteroid tour
trajectory. In particular, the data that must be provided is an ordered
list of the asteroids to be visited. These asteroids must be selected from
the ones available in the problem data for GTOC4. The list of these
asteroids can be found in the webpage of this GTOC edition
(https://sophia.estec.esa.int/gtoc_portal/?page_id=23). Moreover, it is
required to provide a launch date as well as a reasonable initial guess for
the flyby dates of all the asteroids. These have to be provided in Modified
Julian Days [MJD]. All of this data must be stored in a same MATLAB
structure named "Data". This structure must contain the following fields:

Data.Asteroid_Seq_ID ---> The list of the asteroids to be visited. The
asteroids should be listed in the order in which they are to be visited.
Moreover, each asteroid must be referenced through the index number
corresponding to its position in the GTOC4 problem data list. Size: 1D
array.

Data.Launch_Date ---> Earth launch date. Size: scalar [MJD].

Data.Flyby_Dates ---> Initial guess for the flyby dates of each of the
                      asteroids. Size: 1D array [MJD].

ALGORITHM CONTROL PARAMETERS:
These are the different variables that can be used to control the behaviour
of the algorithm. All of these parameters can be changed in the "Algorithm
Options" section at the beginning of the script:

solution_index ---> Flag to select the solution under study. Size: scalar.
                    Five example solutions are already implemented:
                    1 - Moscow
                    2 - Barbee
                    3 - Johnson
                    4 - Original Moscow
                    5 - Impulsive Johson
                    If the user wants to study a self-defined trajectory,
                    this parameter should be set to 6 and the adress of the
                    structure containing the input data must be defined in
                    the "userInitGuess" parameter.

userInitGuess ---> Adress of the file containing input MATLAB structure
                   containing the data of initial guess defined by the
                   user. This parameter is used only if the solution_index
                   parameter has been set to 6. The precise structure that
                   this file must have is the one described in the INPUT
                   DATA part of this description. Type: string.

n_i ---> Number of segments per leg. Size: scalar. (recommended value: 10)

n_simul ---> Size of the optimization windows. Number of arcs that are
             simultaneously optimized. Size: scalar. (recommended value: 2)

asteroidsToRemove ---> Vector containing the index of the asteroids that
                       are to be removed from the trajectory initial guess.
                       Size: 1D array.

kmax ---> Maximum number of MBH local runs (recommended value: 50).

max_rep ---> Maximum number of MBH local runs without improvement
             (recommended value: 5).

rho ---> Maximum percentage of change in MBH (recommended value: 0.15).

resultsAdress ---> Folder adress in which the output file is to be stored.
                   Type: string.

resultsFileName ---> Name of the output file. Type: string.


OUTPUT:
The output of the algorithm is the structure called "sims_flanagan_
trajectory" which contains the impulses to be provided to the spacecraft,
the times in which this impulses are provided and some extra information
needed to reconstruct the solution. This structure contains the following
fields:

Arc# ---> As many fields of this type as legs or arcs are in the
          trajectory. Each of the fields contains two subfields:
          Arc#.vector ---> The cartesian components (x,y,z) of the
                           different impulses to be applied to the S/C in
                           the leg. Size: n_i x 3 [km/s].
          Arc#.times ---> The times in which the impulses are applied.
                          Each entry correspond to the impulse on 
                          Arc#.vector with same index number and represents
                          the elapsed seconds since the beginning of the
                          arc. Size: n_i x 1 [s].

dateModifiers ---> The number of days each flyby has been shifted due to
                   the trajectory optimization, compared to the initial
                   guess. To obtain the optimized flyby dates, it suffices
                   with adding this vector to the provided initial guess
                   for the flyby dates (solution.Flyby_Dates + sims_
                   flanagan_trajectory.dateModifiers). Size: 1 x n_arcs
                   [days].

vInf ---> Cartesian components of the earth departure excess velocity. To
          obtain the intial velocity of the S/C wrt the Sun, this velocity
          should be added to the Earth velocity at the launch date. Size:
          1 x 3 [km/s].

removedAsteroids ---> The index of the asteroids that have been removed
                      from the provided initial guess asteroid sequence. It
                      is equal to the control parameter
                      "asteroidsToRemove". Size: 1D array.

optTime ---> The time taken by the algorithm to obtain the solution. Size:
             scalar [s].

%}

clear; close all; clc;

%% Algorithm Options

%Select solution to be converted (set to 6 for a new trajectory defined by
%the user)
solution_index = 6;

%Adress of the file containing the user's initial guess data (only used if
%"solution_index" has been set to 6)
userInitGuess = "../1.GTOC_Data/TrajectoryInitGuess/EDep12_Tour_10.mat";

%Number of segments per arc in Sims-Flanagan transcription
n_i = 10;

%Size of the optimization window (number of arcs)
n_simul = 2;

%Asteroids to remove
asteroidsToRemove = [];

%MBH parameters
kmax = 50; %Maximum number of local runs
max_rep = 5; %Maximum number of local runs without improvement
rho = 0.15; %Maximum percentage of change

%Adress to store the results and file name
resultsAdress = "../3.Results/Original_Moscow";
resultsFileName = "Original_Moscow.mat";

%% Load GTOC data and trajectory initial guess

addpath("../0.Auxiliar")
addpath("../0.Auxiliar/Basic_2_Body")
addpath("../0.Auxiliar/Transcription_Specifics")
addpath("../0.Auxiliar/Transformations")
addpath("../0.Auxiliar/1-ISF_Transcription")
addpath("../0.Auxiliar/SF_Multiple_Impulses")
addpath("../0.Auxiliar/MBH")

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
solution = removeAsteroids(solution, asteroidsToRemove);

%Loading txt file with GTOC4 Asteroid Data
GTOC4Asteroids = readtable("../1.GTOC_Data/GTOC4.txt");
Asteroids_OE = GTOC4Asteroids(solution.Asteroid_Seq_ID,3:8);

%Earth orbital elements on reference date
pars.OE_Earth_ref = [0.999988049532578, 0.01671681163160, 0.0008854353079654,...
                    175.40647696473, 287.61577546182, 257.60683707535];

%Maximum excess velocity from launch
pars.maxVInf = 4; %[km/s]

%Initial and final dates of the launch window
pars.InitLW = 57023; %[MJD]
pars.FinalLW = 61041; %[MJD]

%Maximum mission time
pars.MaxMissionTime = 10*365; %[MJD]

%% Earth Launch

%Earth position and velocity on the launch date
[r_earth,v_earth] = CartesianState(pars.OE_Earth_ref, (solution.Launch_Date-pars.Epoch.Earth_Ephem)*86400, pars);

%First asteroid orbital elements on reference date (54800 MJD)
Asteroid1_OE = table2array(Asteroids_OE(1,:));

%Position of the first asteroid on S/C flyby date
r_asteroid = CartesianState(Asteroid1_OE, (solution.Flyby_Dates(1)-pars.Epoch.Asteroid_Ephem)*86400, pars);

%Variable to store the ToF of the arcs
ToFs = zeros(length(solution.Flyby_Dates),1);

%Time of flight of the arc earth-first asteroid (in seconds)
ToFs(1) = (solution.Flyby_Dates(1) - solution.Launch_Date)*86400;

%% Asteroid flybys position vectors and TOFs
%{
These are all the positions and ToFs obtained with the initial guess
provided for the flyby dates. As these flyby dates will later be optimized,
the values computed here might be different in the final solution.
%}

%Vector to store the asteroids position vectors
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
for arc_i=2:length(ToFs)
    ToFs(arc_i) = dates(arc_i)-dates(arc_i-1);
end

%Velocity of the final asteroid (used to meet the rendez-vous condition)
Asteroid_OE = table2array(Asteroids_OE(n_asteroids,:));
date = (solution.Flyby_Dates(n_asteroids)-pars.Epoch.Asteroid_Ephem)*86400;
[~,v_finalAsteroid] = CartesianState(Asteroid_OE, date, pars);

%% Sims-Flanagan Transcription Optimization

%Begin chrono
tic;

%Number of arcs in the solution
n_arcs = length(solution.Flyby_Dates);

%Size of the dates optimization window
nArcsDateOpt = n_arcs;

%Mutiple arcs trajectory initial position and mass
r_vec = r_earth;
SCMass = pars.SC.m0;

%Variable to store the impulses
impulses = zeros(n_i,3,n_arcs);

%Variable to store the state vector at each flyby
rv_history = zeros(n_arcs,6);

%Variable to store the masss at each arc
arc_mass = zeros(1,n_arcs);
arc_mass(1) = SCMass;

%Structure to store the modified dates solution
modSolution = solution;

%Variable to store the modified times of flight
Upd_ToFs = ToFs;

%Variable to store the modified flyby positions
Upd_r_flybys = r_flybys;

%Optimize the dates of the full trajectory and launch velocity using the
%1-ISF transcription
modSolution = One_ISF_ModifyDatesWithLaunch(n_arcs,Asteroids_OE,modSolution,pars);

%Optimized launch velocity
v_vec = v_earth + modSolution.vInf;

%Update the times of flight vector with the new flyby dates
Upd_ToFs(1) = (modSolution.Flyby_Dates(1)-modSolution.Launch_Date)*86400;
for arc_k=2:n_arcs
    Upd_ToFs(arc_k) = (modSolution.Flyby_Dates(arc_k)-modSolution.Flyby_Dates(arc_k-1))*86400;
end

%Update the flyby positions for the new flyby dates
for flyby=1:n_arcs
    Upd_r_flybys(flyby,:) = CartesianState(table2array(Asteroids_OE(flyby,:)), (modSolution.Flyby_Dates(flyby)-pars.Epoch.Asteroid_Ephem)*86400, pars);
end

%Optimize the first n_simul-1 arcs in an isolated way to obtain initial
%guess for the first optimization window.
for arc_i=1:n_simul-1
    % Check if the dates optimization window exceeds the legth of the
    % solution
    if nArcsDateOpt > n_arcs-arc_i+1
        actual_nArcsOpt = n_arcs-arc_i+1;
    else
        actual_nArcsOpt = nArcsDateOpt;
    end

    %Optimize the dates of the remaining arcs using the 1-ISF transcription
    modSolution = One_ISF_ModifyDates(arc_i,v_vec,arc_mass(arc_i),actual_nArcsOpt,Asteroids_OE,modSolution,pars);

    %Update the times of flight
    for arc_k=arc_i:arc_i+actual_nArcsOpt-1
        if arc_k == 1
            Upd_ToFs(arc_k) = (modSolution.Flyby_Dates(arc_k)-modSolution.Launch_Date)*86400;
        else
            Upd_ToFs(arc_k) = (modSolution.Flyby_Dates(arc_k)-modSolution.Flyby_Dates(arc_k-1))*86400;
        end
    end

    %Update the flyby positions
    for flyby=arc_i:arc_i+actual_nArcsOpt-1
        Upd_r_flybys(flyby,:) = CartesianState(table2array(Asteroids_OE(flyby,:)), (modSolution.Flyby_Dates(flyby)-pars.Epoch.Asteroid_Ephem)*86400, pars);
    end

    %Save position and velocity of the S/C at departure asteroid
    rv_history(arc_i,:) = [r_vec,v_vec];

    %Time of flight of the arc
    ToF = Upd_ToFs(arc_i);

    %Get the position of the arrival asteroid at the flyby date
    r2_vector = Upd_r_flybys(arc_i,:);

    %Obtain the maximum impulse per segment
    maxdV = pars.SC.T*(ToF/n_i)/(arc_mass(arc_i)*1000);

    %Generate initial guess for the current arc (based on lambert impulses)
    init_guess = SingleArc_InitGuess(r_vec,v_vec,r2_vector,ToF,...
                                n_i,pars).vector;

    %Compute optimal final velocity of the arc using 1-ISF transcription
    VI_post = One_ISF_obtainDesiredFV(arc_i,r_vec,v_vec,ToF,r2_vector,arc_mass(arc_i),Asteroids_OE,modSolution,pars);

    %Optimize SF transcription of the arc to minimize the difference
    %between its final velocity and the computed optimal one
    %(SQP algorithm)
    sims_flanagan_arc = FLI_Setup_NLP(r_vec,v_vec,r2_vector,VI_post,ToF,n_i,pars,init_guess,maxdV);

    %Store the impulses of the SF transcription
    impulses(:,:,arc_i) = sims_flanagan_arc.vector;

    %Update the mass based on the provided impulses
    arc_mass(arc_i+1) = updateArcMass(sims_flanagan_arc.vector.*1000,arc_mass(arc_i),pars);

    %Check that the mass is not lower than the S/C empty mass
    %(otherwise assign it this empty mass value)
    if arc_mass(arc_i+1) < pars.SC.m0 - pars.SC.mprop
        arc_mass(arc_i+1) = pars.SC.m0 - pars.SC.mprop;
    end

    %Propagate the arc to obtain initial position and velocity of next arc
    [r_vec, v_vec] = propagate_segmented(r_vec,v_vec,sims_flanagan_arc.vector,sims_flanagan_arc.times,ToF,pars);
end

%Optimization loop of the optimization windows
for arc_i=n_simul:n_arcs-1
    %Optimize the dates of the remaining arcs using the 1-ISF transcription
    if nArcsDateOpt > n_arcs-arc_i+1
        actual_nArcsOpt = n_arcs-arc_i+1;
    else
        actual_nArcsOpt = nArcsDateOpt;
    end
    modSolution = One_ISF_ModifyDates(arc_i,v_vec,arc_mass(arc_i),actual_nArcsOpt,Asteroids_OE,modSolution,pars);

    %Update the times of flight
    for arc_k=arc_i:arc_i+actual_nArcsOpt-1
        if arc_k == 1
            Upd_ToFs(arc_k) = (modSolution.Flyby_Dates(arc_k)-modSolution.Launch_Date)*86400;
        else
            Upd_ToFs(arc_k) = (modSolution.Flyby_Dates(arc_k)-modSolution.Flyby_Dates(arc_k-1))*86400;
        end
    end

    %Update the flyby positions
    for flyby=arc_i:arc_i+actual_nArcsOpt-1
        Upd_r_flybys(flyby,:) = CartesianState(table2array(Asteroids_OE(flyby,:)), (modSolution.Flyby_Dates(flyby)-pars.Epoch.Asteroid_Ephem)*86400, pars);
    end

    %Save position and velocity of the S/C at departure asteroid
    rv_history(arc_i,:) = [r_vec,v_vec];

    %Time of flight of the arc
    ToF = Upd_ToFs(arc_i);

    %Get the position of the arrival asteroid at the flyby date
    r2_vector = Upd_r_flybys(arc_i,:);

    %Obtain the maximum impulse per segment
    maxdV = pars.SC.T*(ToF/n_i)/(arc_mass(arc_i)*1000);

    %Generate initial guess for the new added arc (based on lambert impulses)
    init_guess = SingleArc_InitGuess(r_vec,v_vec,r2_vector,ToF,...
                                n_i,pars).vector;

    %Compute optimal final velocity of the optimization window using 1-ISF transcription
    VI_post = One_ISF_obtainDesiredFV(arc_i,r_vec,v_vec,ToF,r2_vector,arc_mass(arc_i),Asteroids_OE,modSolution,pars);

    %Optimize the last arc ot the current optimization window using the
    %single arcs approach to use it as initial guess for simultaneous optimization
    %(SQP algorithm)
    sims_flanagan_arc = FLI_Setup_NLP(r_vec,v_vec,r2_vector,VI_post,ToF,n_i,pars,init_guess,maxdV);

    %Store the impulses of the SF transcription
    impulses(:,:,arc_i) = sims_flanagan_arc.vector;

    %Compute the parameters of the optimization window
    r_ini = rv_history(arc_i-(n_simul-1),1:3); %Initial position of the optimization window
    v_ini = rv_history(arc_i-(n_simul-1),4:6); %Initial velocity of the optimization window
    r_obj = Upd_r_flybys(arc_i+1-n_simul:arc_i,:); %Extract the position of the flybys in the optimization window
    arc_ToFs = Upd_ToFs(arc_i+1-n_simul:arc_i); %Extract the time of flight of the arcs in the optimization window
    initMass = arc_mass(arc_i+1-n_simul); %Extract the SC mass at the beginning of the optimization window
    comb_guess = impulses(:,:,arc_i-(n_simul-1):arc_i); %Initial guess for the impulses

    %Optimization of the optimization window (MBH algorithm)
    combined_sims_flanagan = Simultaneous_Setup_MBH(r_ini,v_ini,r_obj,VI_post,arc_ToFs,initMass,comb_guess,n_i,pars.SC.T,pars,kmax,rho,max_rep);

    %Store the results
    impulses(:,:,arc_i-(n_simul-1):arc_i) = combined_sims_flanagan.vector;
    times(:,arc_i-(n_simul-1):arc_i) = combined_sims_flanagan.times;

    %Save the results obtained for the initial arc of the current optimization
    %window (these results are freezed now)
    sims_flanagan_trajectory.(sprintf("Arc%d",arc_i-n_simul+1)).vector = impulses(:,:,arc_i-n_simul+1);
    sims_flanagan_trajectory.(sprintf("Arc%d",arc_i-n_simul+1)).times = times(:,arc_i-n_simul+1);

    %Propagate the arcs to update the state vectors
    rv_updated = propagate_updated(n_simul,r_ini,v_ini,combined_sims_flanagan.vector,combined_sims_flanagan.times,arc_ToFs,pars);
    for j=1:n_simul-1
        rv_history(arc_i-n_simul+1+j,:) = rv_updated(j,:);
    end
    r_vec = rv_updated(end,1:3);
    v_vec = rv_updated(end,4:6);

    %Update the masses according to the provided impulses
    for k=2:n_simul+1
        arc_mass(arc_i-n_simul+k) = updateArcMass(combined_sims_flanagan.vector(:,:,k-1).*1000,arc_mass(arc_i-n_simul+k-1),pars);
        if arc_mass(arc_i-n_simul+k) < pars.SC.m0 - pars.SC.mprop
            arc_mass(arc_i-n_simul+k) = pars.SC.m0 - pars.SC.mprop;
        end
    end

    %Display progress message
    fprintf("Progress. Arc %d finished. \n", arc_i);
end

%Last optimization window (rendez-vous condition: the desired final
%velocity is the one of the final asteroid)
ToF = Upd_ToFs(n_arcs);
r2_vector = Upd_r_flybys(n_arcs,:);
maxdV = pars.SC.T*(ToF/n_i)/(arc_mass(n_arcs)*1000);
init_guess = SingleArc_InitGuess(r_vec,v_vec,r2_vector,ToF,...
                                n_i,pars).vector;
sims_flanagan_arc = FLI_Setup_NLP(r_vec,v_vec,r2_vector,v_finalAsteroid,ToF,n_i,pars,init_guess,maxdV);

impulses(:,:,n_arcs) = sims_flanagan_arc.vector;
r_ini = rv_history(n_arcs-(n_simul-1),1:3); %Initial position of the optimization window
v_ini = rv_history(n_arcs-(n_simul-1),4:6); %Initial velocity of the optimization window
r_obj = Upd_r_flybys(n_arcs+1-n_simul:n_arcs,:); %Extract the position of the flybys in the optimization window
arc_ToFs = Upd_ToFs(n_arcs+1-n_simul:n_arcs); %Extract the time of flight of the arcs in the optimization window
initMass = arc_mass(n_arcs+1-n_simul); %Extract the SC masses at each arc in the optimization window
comb_guess = impulses(:,:,n_arcs-(n_simul-1):n_arcs); %Initial guess for the impulses

combined_sims_flanagan = Simultaneous_Setup_MBH(r_ini,v_ini,r_obj,v_finalAsteroid,arc_ToFs,initMass,comb_guess,n_i,pars.SC.T,pars,kmax,rho,max_rep);

impulses(:,:,n_arcs-(n_simul-1):n_arcs) = combined_sims_flanagan.vector;
times(:,n_arcs-(n_simul-1):n_arcs) = combined_sims_flanagan.times;

%End of the final optimization window
fprintf("Progress. Arc %d finished. \n", n_arcs);

%Save the results in the results structure
for arc_i=n_arcs-n_simul+1:n_arcs
    sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).vector = impulses(:,:,arc_i);
    sims_flanagan_trajectory.(sprintf("Arc%d",arc_i)).times = times(:,arc_i);
end

%Store the date modifiers
sims_flanagan_trajectory.dateModifiers = modSolution.Flyby_Dates - solution.Flyby_Dates;

%Store the obtained v infinity
sims_flanagan_trajectory.vInf = modSolution.vInf;

%Store the removed asteroids
sims_flanagan_trajectory.removedAsteroids = asteroidsToRemove;

%Stop chrono and store the optimization time
executionTime = toc;
sims_flanagan_trajectory.optTime = executionTime;

%% Save the results
SaveAdress = fullfile(resultsAdress, resultsFileName);
save(SaveAdress, "sims_flanagan_trajectory");

%% FUNCTIONS

function rv_history = propagate_updated(n_arcs,r_ini,v_ini,comb_impulses,times,ToFs,pars)

    %Variable to store the state vector
    rv_history = zeros(n_arcs,6);

    %Initial position and velocity
    r_vec = r_ini;
    v_vec = v_ini;
    
    %Loop on the arcs
    for l=1:n_arcs
        %Propagate current arc
        [r_vec, v_vec] = propagate_segmented(r_vec,v_vec,comb_impulses(:,:,l),times(:,l),ToFs(l),pars);

        %Store final position
        rv_history(l,:) = [r_vec, v_vec];
    end

end
