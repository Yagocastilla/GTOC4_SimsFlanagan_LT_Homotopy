function newMass = updateArcMass(impulses, m0, pars)

%Number of impulses
nI = size(impulses,1);

%Compute the total deltaV of the arc
deltaV = 0;
for j = 1:nI
    deltaV = deltaV + norm(impulses(j,:));
end

%Update the mass
newMass = m0*exp(-deltaV/(pars.g0*pars.SC.Isp));

end