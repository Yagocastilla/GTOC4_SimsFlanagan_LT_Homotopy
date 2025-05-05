function E = Kepler(e, M)

%{
This function solves Kepler's equation for the Eccentric Anomaly by means
of a Newton-Raphson iterative algorithm. The algorithm is prepared to
handle only elliptic orbits.
%}

% Initialization
i    = 0;
Nmax = 100;       % maximum number of iterations
Tol  = 1e-10;     % Tolerance
    
% Eccentric Anomaly Initial Guess
E = M;
    
% Newton-Rapshon
while(i < Nmax)
        
K  = E - e*sin(E) - M;
dK = 1 - e*cos(E);
        
    if (abs(K/dK) < Tol)
        break;
    else
        E = E - K/dK;
    end
        
    i = i + 1;

    if i == Nmax
            warning("Newton not converged.")
    end
end

end
    
    
    