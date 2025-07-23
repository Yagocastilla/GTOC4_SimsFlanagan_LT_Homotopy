function [impulses, opt_results] = Simultaneous_Setup_MBH(r_ini,v_ini,...
      r_obj,v_fin,ToFs,initMass,init_guess,n_i,T_limit,pars,kmax,rho,max_rep)

%{
###########################################################################
            Simultaneous Optimization of Several Arcs Setup + MBH
                        Author: Yago Castilla Lamas
###########################################################################

This function takes several arcs (optimization window) of the trajectory
and optimizes them to minimize the maximum mean thrust of the whole
optimization window. The number of arcs in the optimization window is
denoted nsimul The final impulse of each arc is computed by solving a
lambert's problem that automatically satisfies the constraints of the
flyby as in the FLI setup. The optimization problem is solved through the
MBH algorithm.

The algorithm requires a initial guess that has to be provided through the
"init guess" parameter. This initial guess is formed by the impulses
applied in the optimization window in the form of a 3 dimensional matrix of
size n_i x 3 x n_arcs, where the first index is the number of the impulse
in its arc, the second is the cartesian component index (x,y,z) and the
third one is the arc where the impulse is located.

For the scaling of the problem, a non-dimensionalization is performed. The
reference velocity is taken as the mean value of the maximum impulse
magnitudes of each of the arcs of the initial guess in the optimization
window. The reference time is defined as the mean value of the ToF of the
segments in the optimization window and the reference mass is the initial
mass of the optimization window.

INPUTS:
r_ini ---> Initial position vector of the optimization window [km]
v_ini ---> Initial velocity vector of the optimization window [km/s]
r_obj ---> Position vectors of the different flybys in the optimization
           window. Size nsimul x 3 [km]
ToFs ---> Time of flight of the arcs in the optimization window [s]
masses ---> Mass of the SC for each of the arcs (constant) [Kg]
init_guess ---> Initial guess for the optimization window
n_i ---> Number of segments per arc
pars ---> Structure containing problem constants and parameters

OUTPUTS:
impulses.vector ---> Matrix (n_i x 3 x n). Contains the impulses resulting
                     from the optimization.
impulses.times ---> Contains the times at wich the impulses are performed.
                    These times are measured taking as t=0 the beginning of
                    each of the arcs.

###########################################################################
Versions Updates:
V2 ---> Modify the impulses to place them at the center of each segment

%}

%% Problem non-dimensionalization

%Number of arcs to be optimized at the same time
n_arcs = length(ToFs);

%Compute the maximum magnitude of the impulses on each arc of the initial guess
magnitudes = zeros(n_i,n_arcs);
for l=1:n_arcs
    for k=1:n_i
        magnitudes(k,l) = norm(init_guess(k,:,l));
    end
end
MaxImpulses = max(magnitudes,[],1);

%Compute the mean maximum magnitude of the impulses
meanMaxImpulse = sum(abs(MaxImpulses))/n_arcs;

%Compute the mean ToF of the segment of the arcs
Mean_ToF_seg = 0;
for l=1:n_arcs
    Mean_ToF_seg = Mean_ToF_seg + ToFs(l)/n_i;
end
Mean_ToF_seg = Mean_ToF_seg/n_arcs;

%Copy problem data to non-dimensionalize
adim_pars = pars;

%Define reference magnitudes
adim_pars.ref_velocity = meanMaxImpulse;
adim_pars.ref_time = Mean_ToF_seg;
adim_pars.ref_distance = adim_pars.ref_velocity*adim_pars.ref_time;
adim_pars.ref_mass = initMass;
adim_pars.ref_force = adim_pars.ref_mass*adim_pars.ref_distance/(adim_pars.ref_time^2);
adim_pars.mu_sun = pars.mu_sun*adim_pars.ref_time^2/adim_pars.ref_distance^3;
adim_pars.g0 = pars.g0*adim_pars.ref_time^2/adim_pars.ref_distance/1000;
adim_pars.SC.Isp = pars.SC.Isp/adim_pars.ref_time;

%Non-dimensional input variables
r_ini = r_ini/adim_pars.ref_distance;
v_ini = v_ini/adim_pars.ref_velocity;
r_obj = r_obj./adim_pars.ref_distance;
v_fin = v_fin/adim_pars.ref_velocity;
ToFs = ToFs/adim_pars.ref_time;
initMass = initMass/adim_pars.ref_mass;
init_guess = init_guess./adim_pars.ref_velocity;
T_limit = T_limit/1000/adim_pars.ref_force;

%% Optimization Problem: several arcs simultaneous optimization

%Number of free optimization impulses per arc
N = n_i-1;

%Reshape into a vector the initial guess
T_0 = zeros(1,n_arcs*N*3);
dv_0_arc = zeros(1,N*3,1);
SCMass = initMass;
for l=1:n_arcs
    for k=1:N
        dv_0_arc((k-1)*3+1:k*3) = init_guess(k,:,l);
    end
    T_0((l-1)*N*3+1:l*N*3) = dv_0_arc.*SCMass./(ToFs(l)/n_i);
    SCMass = updateArcMass(init_guess(:,:,l),SCMass,adim_pars);
end

%Set objective function
finalVelError = @(T_vector)sim_PostInitError(r_ini,v_ini,r_obj,v_fin,ToFs,initMass,T_vector,adim_pars);

%Set the non-linear inequality constraint
constraint = @(T_vector)thrustConstraint(r_ini,v_ini,r_obj,ToFs,initMass,T_limit,T_vector,adim_pars);

%Penalized objective function
penalObj = @(opt_T_vector)finalVelError(opt_T_vector) + 1e3*max(0, constraint(opt_T_vector));

%Compute the maximum thrust for the initial guess
thrustMagnitudes = zeros(n_i,n_arcs);
SCMass = initMass;
for l=1:n_arcs
    for k=1:n_i
        thrustMagnitudes(k,l) = norm(init_guess(k,:,l))*SCMass/(ToFs(l)/n_i);
    end
    SCMass = updateArcMass(init_guess(:,:,l),SCMass,adim_pars);
end
maxThrust = max(thrustMagnitudes,[],"all");

%Set bounds for the optimization variables
if maxThrust > T_limit
    bounds = maxThrust*ones(1,n_arcs*N*3);
else
    bounds = T_limit*ones(1,n_arcs*N*3);
end

%If not feasible, try a first optimization
if maxThrust > T_limit
    opt_options = optimoptions("fmincon", Algorithm = 'sqp', MaxFunctionEvaluations = 1e5,...
                           MaxIterations = 1e6);
    T_0_x = fmincon(@(x)0,T_0,[],[],[],[],-bounds,bounds,constraint,opt_options);
else
    T_0_x = T_0;
end

%If still not feasible, try with GA
ineq = constraint(T_0_x);
if ineq>1e-3
    GA_finalVelError = @(T_vector)GA_PostInitError(r_ini,v_ini,r_obj,v_fin,ToFs,initMass,T_limit,T_vector,adim_pars);
    ga_options = optimoptions("ga",'Display', 'iter',PlotFcn = @gaplotbestf,InitialPopulationMatrix = T_0_x);
    opt_T_0 = ga(GA_finalVelError,length(T_0),[],[],[],[],-bounds,bounds,[],ga_options);

    %If still not feasible, aditional optimization step
    ineq = constraint(opt_T_0);
    if ineq>1e-3
        opt_options = optimoptions("fmincon", Algorithm = 'sqp', MaxFunctionEvaluations = 1e5,...
                               MaxIterations = 1e6);
        opt_T_0 = fmincon(@(x)0,opt_T_0,[],[],[],[],-bounds,bounds,constraint,opt_options);
    end
else
    opt_T_0 = T_0_x;
end

%Set the local optimization algorithm options
opt_options = optimoptions("fmincon", Algorithm = 'sqp', MaxFunctionEvaluations = 5e4,...
                           MaxIterations = 1e6);

%First local optimization
[opt_solution, fout] = fmincon(penalObj,opt_T_0,[],[],[],[],-bounds,bounds,constraint,opt_options);

%Variable to store current best
opt_results.best_ev = zeros(1,kmax+1);
opt_results.best_constr = zeros(1,kmax+1);

%Monotonic basin hopping loop
j = 0;
rep = 0;
current_best = opt_solution;
current_best_value = fout;
current_best_constraint = constraint(current_best);
opt_results.best_ev(1) = current_best_value;
opt_results.best_constr(1) = current_best_constraint;
while j<kmax
    %Generate the random perturbation vector [-rho, rho]
    pert = zeros(size(opt_solution));
    for l=1:length(pert)
        if current_best(l) == 0
            pert(l) = (2*rand-1)*rho;
        else
            pert(l) = (2*rand-1)*rho*current_best(l);
        end
    end

    %Generate new initial guess
    new_guess = current_best + pert;

    %Local Optimization Algorithm
    [opt_solution, fout] = fmincon(penalObj,new_guess,[],[],[],[],-bounds,bounds,constraint,opt_options);

    %Update counter of local runs
    j = j + 1;

    %If a better solution is found, update the result
    if current_best_constraint<=1e-3
        if fout < current_best_value && constraint(opt_solution)<=1e-3
            current_best = opt_solution;
            current_best_value = fout;
            current_best_constraint = constraint(opt_solution);
            rep = 0;
        else
            rep = rep + 1;
        end
    else
        if constraint(opt_solution) < current_best_constraint
            current_best = opt_solution;
            current_best_value = fout;
            current_best_constraint = constraint(opt_solution);
            rep = 0;
        else
            rep = rep + 1;
        end
    end

    %Save current best
    opt_results.best_ev(j+1) = current_best_value;
    opt_results.best_constr(j+1) = current_best_constraint;

    %After the specified runs, if no improvement is found, stop the algorithm
    if rep == max_rep
        break
    end
end

%Store the results of the optimization
opt_results.local_runs = j+1;
opt_results.repeats = rep;

%% Results

%Variable to store the impulses
impulses.vector = zeros(n_i,3,n_arcs);

%SC mass variable
SCMass = initMass;

%Compute the last impulses of each arc to satisfy lambert's problems
r_vec = r_ini; v_vec = v_ini;
for l=1:n_arcs
    %Reshape the thrusts and convert to impulses
    thrusts_arc = current_best(((l-1)*N*3+1:l*N*3));
    for k=1:N
        impulses.vector(k,:,l) = thrusts_arc((k-1)*3+1:k*3).*(ToFs(l)/n_i)./SCMass;
    end

    %Compute the last impulse of the arc
    impulses.vector(n_i,:,l) = n_last_dV(n_i,r_vec,r_obj(l,:),ToFs(l),v_vec,impulses.vector(1:N,:,l),adim_pars);

    %Time of flight of the segments of the arc
    ToF_segment = ToFs(l)/n_i;

    %Update initial position and velocity
    times = linspace(ToF_segment/2,ToFs(l)-(ToF_segment/2),n_i);
    [r_vec, v_vec] = propagate_segmented(r_vec,v_vec,impulses.vector(1:n_i,:,l),times',ToFs(l),adim_pars);

    %Update the mass
    SCMass = updateArcMass(impulses.vector(:,:,l),SCMass,adim_pars);
end

%Redimensionalize the impulses
for l=1:n_arcs
    for k=1:n_i
        impulses.vector(k,:,l) = impulses.vector(k,:,l)*adim_pars.ref_velocity;
    end
end

%Create a vector with the times of the impulses
impulses.times = zeros(n_i,n_arcs);
for l=1:n_arcs
    for k=1:n_i
        impulses.times(k,l) = (2*k-1)*ToFs(l)*adim_pars.ref_time/(2*n_i);
    end
end

end


%% FUNCTIONS

function finalVel_Diff = sim_PostInitError(r_vec,v_vec,r_obj,v_fin,ToFs,initMass,T_vector,pars)

%Number of arcs, optimization impulses per arc and impulses
n_arcs = length(ToFs); N = length(T_vector)/n_arcs/3; n_i = N + 1;

%Variable to store the impulses
impulses = zeros(n_i,3,n_arcs);

%SC mass variable
SCMass = initMass;

%Value in case of error during execution
errorValue = 10e6;

try
    %Compute the last impulses of each arc to satisfy lambert's problems
    for l=1:n_arcs

        %Reshape the thrusts and convert to impulses
        thrusts_arc = T_vector(((l-1)*N*3+1:l*N*3));
        for k=1:N
            impulses(k,:,l) = thrusts_arc((k-1)*3+1:k*3).*(ToFs(l)/n_i)./SCMass;
        end

        %Time step
        delta_t = ToFs(l)/n_i;
    
        %Loop on the segments
        for k=1:N
            %Propagate the orbit for the first half of the segment
            [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);

            %Compute velocity after the impulse
            v_vec = v_vec + impulses(k,:,l);
        
            %Check if the orbit is hyperbolic
            [~, e, ~, ~, ~, ~] = Build_OE(r_vec, v_vec, pars.mu_sun);
            if e>1
                finalVel_Diff = errorValue;
                return
            end
        
            %Propagate the orbit for the second half of the segment
            [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);
        end

        %Propagate the first half of the last segment of the arc
        [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);
    
        %Compute the last impulse of the arc
        [~,~,~,~,VI,~,~,~] = lambertMR(r_vec, r_obj(l,:), delta_t/2, pars.mu_sun, 0, 0, 0, 1);
        impulses(n_i,:,l) = VI-v_vec;
    
        %Provide last impuse
        v_vec = v_vec + impulses(n_i,:,l);

        %Check if the orbit is hyperbolic
        [~, e, ~, ~, ~, ~] = Build_OE(r_vec, v_vec, pars.mu_sun);
        if e>1
            finalVel_Diff = errorValue;
            return
        end

        %Propagate the second half of the last segment of the arc
        [r_vec, v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);

        %Update the mass
        SCMass = updateArcMass(impulses(:,:,l),SCMass,pars);
    end
    
    %Compute the norm of the difference between the actual final velocity
    %and the desired one
    finalVel_Diff = norm(v_vec - v_fin);

catch
    finalVel_Diff = errorValue;
end

end

%#########################################################################%

function [c,ceq] = thrustConstraint(r_vec,v_vec,r_obj,ToFs,initMass,T_limit,T_vector,pars)

%No equality constraint
ceq = 0;

%Number of arcs, optimization impulses per arc and impulses
n_arcs = length(ToFs); N = length(T_vector)/n_arcs/3; n_i = N + 1;

%Variable to store the impulses
impulses = zeros(n_i,3,n_arcs);

%SC mass variable
SCMass = initMass;

%Value in case of error during execution
errorValue = 10e6;

try
    %Compute the last impulses of each arc to satisfy lambert's problems
    for l=1:n_arcs

        %Reshape the thrusts and convert to impulses
        thrusts_arc = T_vector(((l-1)*N*3+1:l*N*3));
        for k=1:N
            impulses(k,:,l) = thrusts_arc((k-1)*3+1:k*3).*(ToFs(l)/n_i)./SCMass;
        end

        %Time step
        delta_t = ToFs(l)/n_i;
    
        %Loop on the segments
        for k=1:N
            %Propagate the orbit for the first half of the segment
            [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);

            %Compute velocity after the impulse
            v_vec = v_vec + impulses(k,:,l);
        
            %Check if the orbit is hyperbolic
            [~, e, ~, ~, ~, ~] = Build_OE(r_vec, v_vec, pars.mu_sun);
            if e>1
                c = errorValue;
                return
            end
        
            %Propagate the orbit for the second half of the segment
            [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);
        end

        %Propagate the first half of the last segment of the arc
        [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);
    
        %Compute the last impulse of the arc
        [~,~,~,~,VI,~,~,~] = lambertMR(r_vec, r_obj(l,:), delta_t/2, pars.mu_sun, 0, 0, 0, 1);
        impulses(n_i,:,l) = VI-v_vec;
    
        %Provide last impuse
        v_vec = v_vec + impulses(n_i,:,l);

        %Check if the orbit is hyperbolic
        [~, e, ~, ~, ~, ~] = Build_OE(r_vec, v_vec, pars.mu_sun);
        if e>1
            c = errorValue;
            return
        end

        %Propagate the second half of the last segment of the arc
        [r_vec, v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);

        %Update the mass
        SCMass = updateArcMass(impulses(:,:,l),SCMass,pars);
    end
    
    %Compute the magnitude of the impulses
    magnitudes = zeros(n_i,n_arcs);
    for l=1:n_arcs
        for k=1:n_i
            magnitudes(k,l) = norm(impulses(k,:,l));
        end
    end

    %Transform the impulses magnitude into mean thrust magnitude
    thrust = zeros(n_i,n_arcs);
    SCMass = initMass;
    for l=1:n_arcs
        ToF_segment = ToFs(l)/n_i;
        for k=1:n_i
            thrust(k,l) = magnitudes(k,l)*SCMass/ToF_segment;
        end
        SCMass = updateArcMass(impulses(:,:,l),SCMass,pars);
    end
    
    %Compute the maximum impulse
    max_thrust = max(thrust,[],"all");

    %Inequality constraint
    c = max_thrust - T_limit;

catch
    c = errorValue;
end

end

%#########################################################################%

function finalVel_Diff = GA_PostInitError(r_vec,v_vec,r_obj,v_fin,ToFs,initMass,T_limit,T_vector,pars)

%Number of arcs, optimization impulses per arc and impulses
n_arcs = length(ToFs); N = length(T_vector)/n_arcs/3; n_i = N + 1;

%Variable to store the impulses
impulses = zeros(n_i,3,n_arcs);

%SC mass variable
SCMass = initMass;

%Value in case of error during execution
errorValue = 10e6;

try
    %Compute the last impulses of each arc to satisfy lambert's problems
    for l=1:n_arcs

        %Reshape the thrusts and convert to impulses
        thrusts_arc = T_vector(((l-1)*N*3+1:l*N*3));
        for k=1:N
            impulses(k,:,l) = thrusts_arc((k-1)*3+1:k*3).*(ToFs(l)/n_i)./SCMass;
        end

        %Time step
        delta_t = ToFs(l)/n_i;
    
        %Loop on the segments
        for k=1:N
            %Propagate the orbit for the first half of the segment
            [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);

            %Compute velocity after the impulse
            v_vec = v_vec + impulses(k,:,l);
        
            %Check if the orbit is hyperbolic
            [~, e, ~, ~, ~, ~] = Build_OE(r_vec, v_vec, pars.mu_sun);
            if e>1
                finalVel_Diff = errorValue;
                return
            end
        
            %Propagate the orbit for the second half of the segment
            [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);
        end

        %Propagate the first half of the last segment of the arc
        [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);
    
        %Compute the last impulse of the arc
        [~,~,~,~,VI,~,~,~] = lambertMR(r_vec, r_obj(l,:), delta_t/2, pars.mu_sun, 0, 0, 0, 1);
        impulses(n_i,:,l) = VI-v_vec;
    
        %Provide last impuse
        v_vec = v_vec + impulses(n_i,:,l);

        %Check if the orbit is hyperbolic
        [~, e, ~, ~, ~, ~] = Build_OE(r_vec, v_vec, pars.mu_sun);
        if e>1
            finalVel_Diff = errorValue;
            return
        end

        %Propagate the second half of the last segment of the arc
        [r_vec, v_vec] = propagate_kepler(r_vec,v_vec,delta_t/2,pars);

        %Update the mass
        SCMass = updateArcMass(impulses(:,:,l),SCMass,pars);
    end
    
    %Compute the norm of the difference between the actual final velocity
    %and the desired one
    finalVel_Diff = norm(v_vec - v_fin);

    %Compute the magnitude of the impulses
    magnitudes = zeros(n_i,n_arcs);
    for l=1:n_arcs
        for k=1:n_i
            magnitudes(k,l) = norm(impulses(k,:,l));
        end
    end

    %Transform the impulses magnitude into mean thrust magnitude
    thrust = zeros(n_i,n_arcs);
    SCMass = initMass;
    for l=1:n_arcs
        ToF_segment = ToFs(l)/n_i;
        for k=1:n_i
            thrust(k,l) = magnitudes(k,l)*SCMass/ToF_segment;
        end
        SCMass = updateArcMass(impulses(:,:,l),SCMass,pars);
    end

    %Compute the maximum impulse
    max_thrust = max(thrust,[],"all");

    %Inequality constraint
    c = max_thrust - T_limit;

    %Check if the constraint is satisfied (apply the barrier function)
    if c>0
        finalVel_Diff = finalVel_Diff + c*(1e3);
    end

catch
    finalVel_Diff = errorValue;
end

end
