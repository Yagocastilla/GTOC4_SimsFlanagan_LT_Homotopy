function mag2angles = impulseToMag2Angles(r_vec,v_vec,impulse)

%Compute the magnitude of the impulse
impulseMag = norm(impulse);

%Compute steering frame
[u0,v0,w0] = computeSteeringFrame(r_vec,v_vec);

%Compute the direction of the impulse in the steering frame
uT = [dot(impulse,u0),dot(impulse,v0),dot(impulse,w0)]/impulseMag;

%Compute the angles alpha and beta
beta = asin(uT(3));
if cos(beta) == 0
    alpha = 0;
elseif uT(2) < 0
    alpha = -acos(uT(1)/cos(beta));
else
    alpha = acos(uT(1)/cos(beta));
end

%Create the output vector
mag2angles = [impulseMag,alpha,beta];

end