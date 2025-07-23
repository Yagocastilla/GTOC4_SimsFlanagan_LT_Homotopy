function [current_best,opt_results] = MBH(objFunction,constraintFunction,initGuess,infBounds,supBounds,kmax,rho,max_rep,varargin)

    % Default values
    default.onlyFeasible = 0;

    % Treat optional parameters
    if nargin > 8
        onlyFeasible = varargin{1};
    else
        onlyFeasible = default.onlyFeasible;
    end

    %Set the local optimization algorithm options
    opt_options = optimoptions("fmincon", Algorithm = 'sqp', MaxFunctionEvaluations = 5e4,...
                           MaxIterations = 1e6);

    %Variable to store the optimization results
    opt_results.best_ev = zeros(1,kmax+2);
    opt_results.best_constr = zeros(1,kmax+2);

    % Store the initial guess
    opt_results.best_ev(1) = objFunction(initGuess);
    opt_results.best_constr(1) = max(constraintFunction(initGuess));
    
    %First local optimization
    opt_solution = fmincon(objFunction,initGuess,[],[],[],[],infBounds,supBounds,constraintFunction,opt_options);
    
    %If a better solution is found, update the result
    [current_best,current_best_value,current_best_constraint] = updateBestMBH(opt_solution,initGuess,objFunction,constraintFunction);

    %Monotonic basin hopping loop
    j = 0;
    rep = 0;
    opt_results.best_ev(2) = current_best_value;
    opt_results.best_constr(2) = current_best_constraint;
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
        opt_solution = fmincon(objFunction,new_guess,[],[],[],[],infBounds,supBounds,constraintFunction,opt_options);
    
        %Update counter of local runs
        j = j + 1;

        % Store the current best to compare
        oldCurrentBest = current_best;
    
        % Update the best one
        [current_best,current_best_value,current_best_constraint] = updateBestMBH(opt_solution,current_best,objFunction,constraintFunction);

        % If no improvement increase the counter of repetitions
        if current_best == oldCurrentBest
            rep = rep + 1;
        else
            rep = 0;
        end

        %If a better solution is found, update the result
        %{
        optResultConstraint = max(constraintFunction(opt_solution));
        if ~onlyFeasible && current_best_constraint<=1e-3
            if fout < current_best_value && optResultConstraint<=1e-3
                current_best = opt_solution;
                current_best_value = fout;
                current_best_constraint = optResultConstraint;
                rep = 0;
            else
                rep = rep + 1;
            end
        elseif onlyFeasible && optResultConstraint<=1e-3
            current_best = opt_solution;
            current_best_value = fout;
            current_best_constraint = optResultConstraint;
            opt_results.best_ev(j+2) = current_best_value;
            opt_results.best_constr(j+2) = current_best_constraint;
            break
        else
            if optResultConstraint < current_best_constraint
                current_best = opt_solution;
                current_best_value = fout;
                current_best_constraint = optResultConstraint;
                rep = 0;
            else
                rep = rep + 1;
            end
        end
        %}
    
        %Save current best
        opt_results.best_ev(j+2) = current_best_value;
        opt_results.best_constr(j+2) = current_best_constraint;

        % If only feasible activated and feasible solution, terminate the
        % algorithm
        if onlyFeasible && current_best_constraint<=1e-3
            break
        end
    
        %After the specified runs, if no improvement is found, stop the algorithm
        if rep == max_rep
            break
        end
    end
    
    %Store the results of the optimization
    opt_results.local_runs = j+1;
    opt_results.repeats = rep;

end