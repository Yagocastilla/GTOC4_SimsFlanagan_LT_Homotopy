function [u,v,w] = computeSteeringFrame(r_vec,v_vec)

u = v_vec/norm(v_vec);
w = cross(r_vec,v_vec); w = w/norm(w);
v = cross(w,u);

end