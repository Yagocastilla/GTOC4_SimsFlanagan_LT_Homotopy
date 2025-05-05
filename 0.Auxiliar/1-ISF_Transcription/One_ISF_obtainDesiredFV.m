function vFDesired = One_ISF_obtainDesiredFV(arc_i,r_ini,v_ini,ToF,r_fin,initMass,asteroidsOE,solution,pars)
    
    %Number of arcs until the end of the trajectory
    nArcs = length(solution.Flyby_Dates)-arc_i;

    %Aproximate the mass lost on the arc and guess for desired vf through central transcription
    [r_vec,v_vec] = propagate_kepler(r_ini,v_ini,ToF/2,pars);
    [~,~,~,~,VI,VF_guess] = lambertMR(r_vec,r_fin,ToF/2,pars.mu_sun, 0, 0, 0, 1);
    centLambImp = norm(VI - v_vec);
    aproxMass = initMass*exp(-norm(centLambImp)*1000/(pars.g0*pars.SC.Isp));
    if aproxMass < pars.SC.m0 - pars.SC.mprop
        aproxMass = pars.SC.m0 - pars.SC.mprop;
    end

    %Set objective function
    objFunc = @(v_ini)trajectoryMaxThrust(arc_i+1,v_ini,aproxMass,nArcs,asteroidsOE,solution,pars);

    %Optimization algorithm
    opt_options = optimoptions("fmincon", MaxFunctionEvaluations = 1e5,...
                                MaxIterations = 1e6);
    vFDesired = fmincon(objFunc,VF_guess,[],[],[],[],[],[],[],opt_options);

end

%% Functions

function maxThrust = trajectoryMaxThrust(init_arc,v_ini,initMass,nArcs,asteroidsOE,solution,pars)

    %Compute the new dates
    flybyDates = zeros(1,nArcs+1);
    if init_arc == 1
        flybyDates(1) = solution.Launch_Date;
        flybyDates(2:end) = solution.Flyby_Dates(init_arc:init_arc+nArcs-1);
    else
        flybyDates = solution.Flyby_Dates(init_arc-1:init_arc+nArcs-1);
    end

    %Compute the position of the flybys with the modified dates
    n_asteroids = length(flybyDates);
    rFlybys = zeros(n_asteroids,3);
    if init_arc == 1
        %Earth orbital elements on reference date (54000 MJD)
        OE_Earth_0 = [0.999988049532578, 0.01671681163160, 0.0008854353079654,...
            175.40647696473, 287.61577546182, 257.60683707535];
        %Earth position on the launch date
        rFlybys(1,:) = CartesianState(OE_Earth_0, (solution.Launch_Date-pars.Epoch.Earth_Ephem)*86400, pars);
    else
        %Asteroid orbital elements on reference date (54800 MJD)
        Asteroid_OE = table2array(asteroidsOE(init_arc - 1,:));
        %Elapsed time since the reference date until the flyby date in seconds
        date = (flybyDates(1)-pars.Epoch.Asteroid_Ephem)*86400;
        %Position of the asteroid on S/C flyby
        rFlybys(1,:) = CartesianState(Asteroid_OE, date, pars);
    end
    for asteroid=2:n_asteroids
        %Asteroid orbital elements on reference date (54800 MJD)
        Asteroid_OE = table2array(asteroidsOE(init_arc + asteroid - 2,:));
    
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

    %Compute the thrust required for each arc with the central lambert
    %transcription
    centlambImpThrustMag = zeros(1,nArcs);
    VF_pre = v_ini;
    SCMass = initMass;

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