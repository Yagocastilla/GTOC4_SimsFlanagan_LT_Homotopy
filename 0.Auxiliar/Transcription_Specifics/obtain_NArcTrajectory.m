function [trajectory,impulses_r] = ...
         obtain_NArcTrajectory(r_ini,v_ini,ToF,impulses,mesh_refiner,pars)

%{
###########################################################################
                         Obtain N-Arc Trajectory
                       Author: Yago Castilla Lamas
###########################################################################

This function takes an arc and the impulses that this arc may contain. Then
the time is discretized in the desired number of points between each
impulse and the position is obtained in each of the points returning the
trajectory followed by the spacecraft:

INPUTS:
r_ini ---> Initial position vector [km]
v_ini ---> Initial velocity vector [km/s]
ToF ---> Time of flight of the arc [s].
impulses ---> Structure containing two informations about the impulses
              (output of the "Single_Arc_..." functions):
              impulses.vector: Matrix comtaining the impulses of the arc.  
                               The rows are the impulses and the columns   
                               the cartesian components. Size Nx3.
              impulses.times: Times relative to the beginning of the       
                              segment where the impulses take place. Size 
                              Nx1.
mesh_refiner ---> Number of points in which each segment in discretized.
pars ---> Structure containing problem constants and parameters.

OUTPUTS:
trajectory ---> Matrix of size (N x mesh_refiner x 3) containing the
                position vectors of the S/C at each of the discretized
                points.
impulses_r ---> Matrix of size Nx3 containing the position vectors of the
                points where the impulses are performed.

###########################################################################
Versions Updates:
V2 ---> Modify the impulses to place them at the center of each segment

%}

%Number of impulses or segments
n_i = length(impulses.times);

%Time of flight of the segments
ToF_segment = ToF/n_i;

%Variable to store the trajectory of the different segments
segments_xyz = zeros(n_i,mesh_refiner,3); %[arc,point,coordinate]

%Variable to store the position of the impulses
impulses_r = zeros(n_i,3);

%Number of points in the first half of each segment
num_points_first = floor(mesh_refiner/2);

%Compute initial position and velocity of the arc
r_vec = r_ini;
v_vec = v_ini;

%Loop on the segments
for k=1:n_i

    %Obtain the trajectory of the first half of the segment
    segments_xyz(k,1:num_points_first,:) = obtain_Segment_Trajectory(r_vec,v_vec,ToF_segment/2,...
                                            pars,num_points_first);

    %Propagate the orbit for the first half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

    %Store the position of the impulse
    impulses_r(k,:) = r_vec;

    %Apply the impulse
    v_vec = v_vec + impulses.vector(k,:);

    %Obtain the trajectory of the second half of the segment
    segments_xyz(k,num_points_first+1:mesh_refiner,:) = obtain_Segment_Trajectory(r_vec,v_vec,ToF_segment/2,...
                                            pars,mesh_refiner-num_points_first);

    %Propagate the orbit for the second half of the segment
    [r_vec,v_vec] = propagate_kepler(r_vec,v_vec,ToF_segment/2,pars);

end

%Reshape the vector to have a unique column with all the points
trajectory = zeros(n_i*mesh_refiner,3);
for k=1:n_i
    trajectory(mesh_refiner*(k-1)+1:mesh_refiner*k,:) = segments_xyz(k,:,:);
end

end

