function cartesianComponents = convert2cartesianRef(r_vec,v_vec,mag2angles)

%Compute the unit vector direction in the steering reference frame
u_steering = [cos(mag2angles(3))*cos(mag2angles(2)),...
              cos(mag2angles(3))*sin(mag2angles(2)),...
              sin(mag2angles(3))];

%Compute steering frame
[u0,v0,w0] = computeSteeringFrame(r_vec,v_vec);

%Build the rotation matrix from the steering frame to the global
R = [u0',v0',w0'];

%Compute the cartesian components
cartesianComponents = mag2angles(1)*(R*u_steering');

end