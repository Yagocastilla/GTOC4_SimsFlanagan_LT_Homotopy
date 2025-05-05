%% GTOC 4 EXISTING SOLUTIONS %%

clear all; close all; clc;

%% INITIALIZATION %%
% Loading txt file with Asteroid Data
GTOC4Asteroids = readtable('GTOC4.txt');

% Convert Asteroid ID's to a string for manipulation
Asteroid_Names = string(GTOC4Asteroids{:,1});
Asteroid_Names = erase(Asteroid_Names,"'");  % Remove ' for string comparison

%% PARAMETERS FROM THE GTOC 4 PROBLEM %%
% Heliocentric J2000 Orbital Elements for the Selected Asteroids [a,e,i,RAAN,w,MA]
Ast_OE_RefEpoch = table2array(GTOC4Asteroids(:,3:8));  %[AU,-,deg,deg,deg,deg]

% Retrieve Epoch at which Asteroid Orbital Elements are given
pars.Epoch.Asteroid_Ephem = table2array(GTOC4Asteroids(1,2));   %[MDJ]

% Heliocentric Orbital Elements for the Earth @ J2000 [a,e,i,RAAN,w,MA]
OE_Earth_0 = [0.999988049532578, 0.01671681163160, 0.0008854353079654,...
    175.40647696473, 287.61577546182, 257.60683707535];

% Epoch at which Earth Orbital Elements are given
pars.Epoch.Earth_Ephem = 54000;   %[MDJ]

% Defining Problem Constants 
pars.mu_sun = 1.32712440018*10^(11);  %[km^3/s^2]
pars.AU     = 1.49597870691*10^(8);   %[km]
pars.g0     = 9.80665;                %[m/s^2] 
pars.Day    = 86400;                  %[s]
pars.Year   = 365.25;                 %[days]
pars.TOF    = 1;                      %[years] Available Mission Time of Flight
pars.JD     = 2400000.5;

% Defining Spacecraft Parameters
pars.SC.DVLaunch = 4;      %[km/s] Maximum hyperbolic Earth departure velocity
pars.SC.m0       = 1500;   %[kg] Initial Wet mass
pars.SC.mprop    = 1000;   %[kg] Propellant mass
pars.SC.Isp      = 3000;   %[s]  Specific Impulse
pars.SC.T        = 0.135;  %[N]  Maximum Engine Thrust
pars.SC.Accel    = pars.SC.T/ pars.SC.m0; %[m/s2] Engine Acceleration
pars.SC.DV       = -pars.SC.Isp*pars.g0*log(1 - pars.SC.mprop/pars.SC.m0 ); %[m/s] Available DV 

%% RETRIEVE DATA FROM SELECTED GTOC 4 SOLUTION %%
Flag = 5;   % 0 = Moscow Solution // 1 = Johnson Solution // 2 = Barbee Solution
[Data] = Data_GTOC4_Solution(Flag, Asteroid_Names);

