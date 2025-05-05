function modSol = One_ISF_ModifyDatesWithLaunch(nArcs,asteroidsOE,solution,pars)

    %% Setup

    %Earth position on the launch date
    [r_earth,v_earth] = CartesianState(pars.OE_Earth_ref, (solution.Launch_Date-pars.Epoch.Earth_Ephem)*86400, pars);

    if isfield(solution, 'vInf')
        VInfCartesian = solution.vInf;
    else
        %First asteroid orbital elements on reference date (54800 MJD)
        Asteroid1_OE = table2array(asteroidsOE(1,:));
    
        %Position of the first asteroid on S/C flyby date
        r_asteroid = CartesianState(Asteroid1_OE, (solution.Flyby_Dates(1)-pars.Epoch.Asteroid_Ephem)*86400, pars);
    
        %Time of flight of the arc earth-first asteroid
        ToF = (solution.Flyby_Dates(1) - solution.Launch_Date)*86400;
    
        %Compute the lambert's arc
        [~,~,~,~,v_launch] = lambertMR(r_earth, r_asteroid, ToF, pars.mu_sun, 0, 0, 0, 1);
    
        % V infinity
        VInfCartesian = v_launch - v_earth;
    end

    %Transform to mag2imp
    VInf = impulseToMag2Angles(r_earth,v_earth,VInfCartesian);

    %% Optimization
    %Copy the solution
    modSol = solution;

    %Initial guess
    initGuess = [VInf(1)/pars.maxVInf; VInf(2)/pi; VInf(3)/(pi/2); zeros(nArcs-1,1)];

    %Set objective function
    objFunc = @(modifiers)trajectoryMaxThrust(nArcs,asteroidsOE,solution,modifiers,pars);

    %Set the linear inequalities contraints matrices
    A = zeros(nArcs,nArcs-1);
    A(1,1) = -1;
    for k = 2:nArcs-1
        A(k,k-1) = 1;
        A(k,k) = -1;
    end
    A(nArcs,nArcs-1) = 1;
    A = [zeros(nArcs,3),A];
    %A(nArcs+1,nArcs) = 1;
    b = zeros(nArcs,1);
    b(1) = solution.Flyby_Dates(1) - solution.Launch_Date;
    for k = 2:nArcs
        b(k) = solution.Flyby_Dates(k) - solution.Flyby_Dates(k-1);
    end
    %b(nArcs+1) = (solution.Flyby_Dates(init_arc+nArcs) - solution.Flyby_Dates(init_arc))/nArcs;

    %Set bounds
    lb = [0; -1; -1; -Inf*ones(nArcs-1,1)];
    ub = [1; 1; 1; Inf*ones(nArcs-1,1)];

    %Set optimization algorithm options
    opt_options = optimoptions("fmincon", MaxFunctionEvaluations = 1e5,...
                                MaxIterations = 1e6);

    %Optimization Algorithm
    modifiers = fmincon(objFunc,initGuess,A,b,[],[],lb,ub,[],opt_options);

    %Split the result in date modifiers and v infinity
    vInf = [modifiers(1)*pars.maxVInf, modifiers(2)*pi, modifiers(3)*pi/2];
    vInfCartesian = convert2cartesianRef(r_earth,v_earth,vInf);
    dateModifiers = modifiers(4:end);

    %Compute the new dates
    modSol.Flyby_Dates(1:nArcs-1) = modSol.Flyby_Dates(1:nArcs-1) + dateModifiers.';

    %Store the initial vInf
    modSol.vInf = vInfCartesian';
end

function maxThrust = trajectoryMaxThrust(nArcs,asteroidsOE,solution,modifiers,pars)

    %Split optimization variables in date modifiers and v infinity
    vInf = [modifiers(1)*pars.maxVInf, modifiers(2)*pi, modifiers(3)*pi/2];
    dateModifiers = modifiers(4:end);

    %Compute the new dates
    flybyDates = zeros(1,nArcs+1);
    flybyDates(1) = solution.Launch_Date;
    flybyDates(2:end) = solution.Flyby_Dates(1:nArcs);
    flybyDates(2:end-1) = flybyDates(2:end-1) + dateModifiers.';

    %Compute the position of the flybys with the modified dates
    n_asteroids = length(flybyDates);
    rFlybys = zeros(n_asteroids,3);

    %Earth position and velocity on the launch date
    [rFlybys(1,:), vEarth] = CartesianState(pars.OE_Earth_ref, (solution.Launch_Date-pars.Epoch.Earth_Ephem)*86400, pars);

    for asteroid=2:n_asteroids
        %Asteroid orbital elements on reference date (54800 MJD)
        Asteroid_OE = table2array(asteroidsOE(asteroid - 1,:));
    
        %Elapsed time since the reference date until the flyby date in seconds
        date = (flybyDates(asteroid)-pars.Epoch.Asteroid_Ephem)*86400;
        
        %Position of the asteroid on S/C flyby
        rFlybys(asteroid,:) = CartesianState(Asteroid_OE, date, pars);
    end

    %Compute the new times of flight
    ToFs = zeros(1,nArcs);
    for k=1:nArcs
        ToFs(k) = (flybyDates(k+1) - flybyDates(k))*86400;
    end

    %Compute initial velocity of the spacecraft
    vInfCartesian = convert2cartesianRef(rFlybys(1,:),vEarth,vInf);
    v_ini = vEarth + vInfCartesian';

    %Compute the thrust required for each arc with the central lambert
    %transcription
    centlambImpThrustMag = zeros(1,nArcs);
    VF_pre = v_ini;
    SCMass = pars.SC.m0;

    %Value in case of error during execution
    errorValue = 10e6;
    
    try
        for arc_i=1:nArcs
            %Propagate the orbit for half an arc
            [r_vec,v_vec] = propagate_kepler(rFlybys(arc_i,:),VF_pre,ToFs(arc_i)/2,pars);
        
            %Compute Lambert's Arc
            [~,~,~,~,VI,VF] = lambertMR(r_vec,rFlybys(arc_i+1,:),ToFs(arc_i)/2,pars.mu_sun, 0, 0, 0, 1);
        
            %Compute Lmabert's Impulse
            centlambImp = VI - v_vec;
        
            %Convert it to thrust
            centlambImpThrustMag(arc_i) = norm(centlambImp)*1000*SCMass/ToFs(arc_i);
        
            %Update velocity
            VF_pre = VF;

            %Update the mass
            SCMass = SCMass*exp(-norm(centlambImp)*1000/(pars.g0*pars.SC.Isp));
            if SCMass < pars.SC.m0 - pars.SC.mprop
                SCMass = pars.SC.m0 - pars.SC.mprop;
            end
        end
    
        %Obtain the maximum thrust
        maxThrust = max(centlambImpThrustMag);
    catch
        maxThrust = errorValue;
    end
end

