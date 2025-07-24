function [halfWayDSMs,halfWayDSMThrusts,halfWayDSMMasses] = computeHalfWayDSM(rFlybys,v_ini,initMass,ToFs,pars)

    % Compute number of legs
    nLegs = length(ToFs);
    
    % Check that the number of targets match
    if nLegs ~= size(rFlybys,1)-1
        error("Sizes of rFlybys and ToFs are not consistent.")
    end
    
    %Compute the thrust required for each arc with the Half-way DSM
    %transcription
    halfWayDSMs = zeros(nLegs,3);
    halfWayDSMThrusts = zeros(nLegs,3);
    halfWayDSMMasses = zeros(1,nLegs);
    VF_pre = v_ini;
    SCMass = initMass;
    
    % Loop in the legs
    for leg_i=1:nLegs
        % Store the mass for the current leg
        halfWayDSMMasses(leg_i) = SCMass;
    
        %Propagate the orbit for half an arc
        [r_vec,v_vec] = propagate_kepler(rFlybys(leg_i,:),VF_pre,ToFs(leg_i)/2,pars);
    
        %Compute Lambert's Arc
        [~,~,~,~,VI,VF] = lambertMR(r_vec,rFlybys(leg_i+1,:),ToFs(leg_i)/2,pars.mu_sun, 0, 0, 0, 1);
    
        %Compute half-way DSM
        halfWayDSMs(leg_i,:) = VI - v_vec;
    
        %Convert it to thrust
        halfWayDSMThrusts(leg_i,:) = halfWayDSMs(leg_i,:)*1000*SCMass/ToFs(leg_i);
    
        %Update velocity
        VF_pre = VF;
    
        %Update the mass
        SCMass = SCMass*exp(-norm(halfWayDSMs(leg_i,:))*1000/(pars.g0*pars.SC.Isp));
        if SCMass < pars.SC.m0 - pars.SC.mprop
            SCMass = pars.SC.m0 - pars.SC.mprop;
        end
    end

end