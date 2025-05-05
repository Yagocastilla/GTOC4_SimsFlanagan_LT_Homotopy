function impulses = FLI_Setup_NLP(r_ini,v_ini,r_fin,v_fin,ToF,n_i,pars,init,dv_limit)

%{
###########################################################################
               Single Arc Final Lambert's Impulse Setup + NLP
                        Author: Yago Castilla Lamas
###########################################################################

This function takes an arc and divides it in N equally spaced segments in
time that still satisfies the lambert's problem. An impulse is placed at
the beginning of each segment. The last segment is not treated as an
optimization variable and instead is computed by solving a lambert's
problem that automatically satisfies the constraint in the final position.
The programme minimizes the maximum impulse of the arc. The function
returns the N impulses as well as the times at which these are performed.

Previously to the optimization, the problem is non-dimensionalized using
the lambert's impulse of the arc divided by the selected number of impulses
as the reference velocity.

An initial guess have to be provided to initialize the optimization
algorithm. This can be done through the "init" parameter.

The optimization problem is solved through the fmincon non-linear
programming algorithm.

INPUTS:
r_ini ---> Initial position vector [km]
v_ini ---> Initial velocity vector [km/s]
r_fin ---> Final position vector [km]
ToF ---> Time of flight of the arc [s]
n_i ---> Number of impulses or segments
pars ---> Structure containing problem constants and parameters
init ---> Initial guess. Matrix (n_i x 3). The rows are the different
          impulses and the columns correspond to the 3 cartesian 
          components.
maxdV ---> Maximum impulse that the propulsion system can provide in a
           segment (only needed when using initial guess 4) [km/s]

OUTPUTS:
impulses.vector ---> Matrix (n_i x 3). The rows are the different impulses
                     and the columns correspond to the 3 cartesian
                     components.
impulses.times ---> Vector containing the times at wich the impulses are
                    performed. These times are measured taking as t=0 the
                    beginning of the arc.

###########################################################################
Versions Updates:
V2 ---> Modify the impulses to place them at the center of each segment

%}

%Check that the number of impulses and the initial guess are consistent
if n_i ~= size(init,1)
    error("The number of impulses of the arc and the initial guess" + ...
        " are not consistent.")
end

%% Problem non-dimensionalization

%Compute lambert's impulse
[~,~,~,~,VI_lamb] = lambertMR(r_ini, r_fin, ToF, pars.mu_sun, 0, 0, 0, 1);
lamb_impulse = VI_lamb - v_ini;

%Copy problem data to non-dimensionalize
adim_pars = pars;

%Define reference magnitudes
adim_pars.ref_distance = norm(r_ini);
adim_pars.ref_velocity = norm(lamb_impulse)/n_i;
if norm(lamb_impulse) < 0.1*dv_limit
    adim_pars.ref_velocity = dv_limit;
end
adim_pars.ref_time = adim_pars.ref_distance/adim_pars.ref_velocity;
adim_pars.ref_mu = (adim_pars.ref_distance^3)/(adim_pars.ref_time^2);

%Non-dimensional input variables
r_ini = r_ini/adim_pars.ref_distance;
v_ini = v_ini/adim_pars.ref_velocity;
r_fin = r_fin/adim_pars.ref_distance;
v_fin = v_fin/adim_pars.ref_velocity;
ToF = ToF/adim_pars.ref_time;
adim_pars.mu_sun = pars.mu_sun/adim_pars.ref_mu;
dv_limit = dv_limit/adim_pars.ref_velocity;

%% Optimization Problem: N Impulses, minimize maximum dV

%Set the objective function
finalVelError = @(opt_dv_vector)PostInitError(n_i,r_ini,r_fin,ToF,v_ini,v_fin,opt_dv_vector,adim_pars);

%Set the non-linear inequality constraint
constraint = @(opt_dv_vector)thrustConstraint(n_i,r_ini,r_fin,ToF,v_ini,opt_dv_vector,adim_pars,dv_limit);

%Penalized objective function
penalObj = @(opt_dv_vector)finalVelError(opt_dv_vector) + 1e3*max(0, constraint(opt_dv_vector));

%Reshape intial guess into a vector
dv_0 = init(1:n_i-1,:);
opt_dv_0 = zeros(1,size(dv_0,1)*3);
for k=1:size(dv_0,1)
    opt_dv_0((k-1)*3+1:k*3) = dv_0(k,:);
end

%Set bounds for the optimization variables
max_init = n_max_dV(n_i,r_ini,r_fin,ToF,v_ini,dv_0,adim_pars);
if max_init>dv_limit
    bounds = ones(size(opt_dv_0))*max_init;
else
    bounds = ones(size(opt_dv_0))*dv_limit;
end

% Genetic algorithm for initial guess if lambert's impulse is too small
%ga_options = optimoptions("ga",'Display', 'iter',PlotFcn = @gaplotbestf);
ga_options = optimoptions("ga");
if norm(lamb_impulse)/adim_pars.ref_velocity < 0.01*dv_limit
    opt_dv_0 = ga(finalVelError,length(opt_dv_0),[],[],[],[],-bounds,bounds,constraint,ga_options);
end

%Set optimization algorithm options
opt_options = optimoptions("fmincon", Algorithm = 'sqp',...
                           MaxFunctionEvaluations = 1e5,...
                           MaxIterations = 1e6);

%If the init guess is not feasible run a first optimization
if constraint(opt_dv_0)>1e-3
    opt_dv_0 = fmincon(@(x)0,opt_dv_0,[],[],[],[],-bounds,bounds,constraint,opt_options);
end

%Optimization Algorithm
opt_solution = fmincon(penalObj,opt_dv_0,[],[],[],[],-bounds,bounds,constraint,opt_options);

%Reshape the result
n_dV_opt = zeros(size(dv_0));
for k=1:size(n_dV_opt,1)
    n_dV_opt(k,:) = opt_solution((k-1)*3+1:k*3);
end

%% Results

%Obtain the required final impulse
last_dV = n_last_dV(n_i,r_ini,r_fin,ToF,v_ini,n_dV_opt,adim_pars);

%Create a vector with the impulses and with original units
impulses.vector = [n_dV_opt;last_dV].*adim_pars.ref_velocity;

%Create a vector with the times of the impulses and with original units
impulses.times = zeros(size(impulses.vector,1),1);
for j=1:length(impulses.times)
    impulses.times(j) = (2*j-1)*ToF*adim_pars.ref_time/(2*n_i);
end

end


%% FUNCTIONS

%#########################################################################%
%{
This function takes the optimization variables (the components of the 
impulses arranged in a one-dimensional vector). It reshape them into a two
dimensional array with a row for each impulse and a column for each
component, and computes the maximum impulse magnitude of the arc using the
"n_max_dV" function. An error can occur during the execution because the
optimization variables take values that result in very extreme orbits and
manoeuvres for which the subroutines of the 2-body problem are not
prepared. In case of such an error, the output is assigned a prohivitive
value.
%}

function finalVel_Diff = PostInitError(n_i,r_ini,r_fin,ToF,v_ini,v_fin,opt_dv_vector,pars)

%Reshape the optimization variables
dv_vector = zeros(n_i,3);
for k=1:n_i-1
    dv_vector(k,:) = opt_dv_vector((k-1)*3+1:k*3);
end

%Value in case of error during execution
errorValue = 10e6;

try
    %Compute the last impulse
    dv_vector(n_i,:) = n_last_dV(n_i,r_ini,r_fin,ToF,v_ini,dv_vector(1:n_i-1,:),pars);

    %Create a vector with the times of the impulses
    times = zeros(n_i,1);
    for j=1:length(times)
        times(j) = (2*j-1)*ToF/(2*n_i);
    end

    %Compute the final velocity
    [~, v_vec] = propagate_segmented(r_ini,v_ini,dv_vector,times,ToF,pars);

    %Compute the norm of the difference between the actual final velocity
    %and the desired one
    finalVel_Diff = norm(v_vec - v_fin);
catch
    finalVel_Diff = errorValue; %In case of error
end

end

%#########################################################################%

%{
This function takes the reshaped optimization impulses. It computes the
required last impulse by solving the lambert's problem and computes the
maximum impulse magnitude of the arc.
%}

function max_dV = n_max_dV(n,r_ini,r_fin,ToF,v_ini,dv_vector,pars)
%Compute the last impulse
last_dV = n_last_dV(n,r_ini,r_fin,ToF,v_ini,dv_vector,pars);

%Compute the magnitude of each impulse
impulse_mag = zeros([1,size(dv_vector,1)+1]);
for k=1:size(dv_vector,1)
    impulse_mag(1,k) = norm(dv_vector(k,:));
end
impulse_mag(1,end) = norm(last_dV);

%Compute the maximum impulse of the arc
max_dV = max(impulse_mag);

end

%#########################################################################%

function [c,ceq] = thrustConstraint(n_i,r_ini,r_fin,ToF,v_ini,opt_dv_vector,pars,dv_limit)

%No equality constraint
ceq = 0;

%Reshape the optimization variables
dv_vector = zeros(n_i-1,3);
for k=1:size(dv_vector,1)
    dv_vector(k,:) = opt_dv_vector((k-1)*3+1:k*3);
end

%Value in case of error during execution
errorValue = 10e6;

try
    %Compute the maximum impulse
    max_dV = n_max_dV(n_i,r_ini,r_fin,ToF,v_ini,dv_vector,pars);

    %Inequality constraint
    c = max_dV - dv_limit;
catch
    c = errorValue;
end
end
