function [r_vector, v_vector] = CartesianState(Body_OE_RefEpoch, Elapsed_Time, pars)

 % This function obtains the cartesian position and velocity vector of an
 % asteroid at a certain Epoch, converting from orbital elements.


% Author: José Carlos Garcia
% Last revision: 06/03/2023 (VERIFIED USING JPL HORIZONS SYSTEM)

%% INPUTS %%
%
% - Body_OE_RefEpoch: Matrix containing the Orbital Elements of the Body at the
%                     Reference Epoch. Dimensions = [1 x 6]. Orbital elements are
%                     expressed in Heliocentric J2000 frame and the columns are
%                     ordered as: [a,e,i,RAAN,w,MA]; with units [AU,-,deg,deg,deg,deg]
%
% - Elapsed_Time: Time difference between the Reference Epoch and the Epoch at which
%                 the state vector has to be obtained. Units: [seconds]
%
% - pars: Structure containing problem constants and parameters

%% OUTPUTS %% 
%
% - r_vector: Position vector in [km] at Epoch in Cartesian coordinates. Dimensions = [1 x 3]
%
% - v_vector: Velocity vector in [km/s] at Epoch in Cartesian coordinates. Dimensions = [1 x 3]


%% OBTAINING MEAN ANOMALY OF THE BODY AT DESIRED EPOCH %%
a       = Body_OE_RefEpoch(1)*pars.AU;  %[km]

% MA_Body = MA_Body_RefEpoch(idx) + sqrt(pars.mu_sun/(a^3))*Elapsed_Time;
MA_Body = deg2rad(Body_OE_RefEpoch(6)) + sqrt(pars.mu_sun/(a^3))*Elapsed_Time;


% Retaining Mean Anomaly inside [0,2pi] range
MA_Body =  (MA_Body/(2*pi) - floor(MA_Body/(2*pi)))*(2*pi);   % [rad]


%% SOLVING KEPLER's EQUATION FOR ECCENTRIC ANOMALY %%
e = Body_OE_RefEpoch(2);

E = Kepler(e, MA_Body); %[rad]

%% SOLVING FOR TRUE ANOMALY %%

TA = 2*atan(sqrt((1+e)/(1-e))* tan(E/2));       %[rad]


%% OBTAINING VELOCITY MAGNITUDE & FLIGHT PATH ANGLE %%
r = (a*(1 - e^2))/(1 + e*cos(TA));                    %[km]

v_norm = sqrt((2*pars.mu_sun/r) - (pars.mu_sun/a));   %[km/s]

gamma  = atan((e*sin(TA))/(1 + e*cos(TA)));           %[rad]


%% COMPUTING STATE VECTOR %%

inc  = deg2rad(Body_OE_RefEpoch(3));    %[rad] Body Inclination
RAAN = deg2rad(Body_OE_RefEpoch(4));    %[rad] Body Longitude of Ascending Node
w    = deg2rad(Body_OE_RefEpoch(5));    %[rad] Body Argument of Periapsis

r_x = r*(cos(TA + w)*cos(RAAN) - sin(TA + w)*cos(inc)*sin(RAAN));    %[km]
r_y = r*(cos(TA + w)*sin(RAAN) + sin(TA + w)*cos(inc)*cos(RAAN));    %[km]
r_z = r*( sin(TA + w)*sin(inc));                                     %[km]

r_vector = [r_x, r_y, r_z];        %[km]

v_x = v_norm*(-sin(TA + w - gamma)*cos(RAAN) - cos(TA + w - gamma)*cos(inc)*sin(RAAN));    %[km/s]
v_y = v_norm*(-sin(TA + w - gamma)*sin(RAAN) + cos(TA + w - gamma)*cos(inc)*cos(RAAN));    %[km/s]
v_z = v_norm*(cos(TA + w - gamma)*sin(inc));                                               %[km/s]

v_vector =  [v_x, v_y, v_z];        %[km/s]

end
