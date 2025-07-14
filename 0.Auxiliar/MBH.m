function [current_best,opt_results] = MBH(objFunction,constraintFunction,initGuess,infBounds,supBounds,kmax,rho,max_rep)

    %Set the local optimization algorithm options
    opt_options = optimoptions("fmincon", Algorithm = 'sqp', MaxFunctionEvaluations = 5e4,...
                           MaxIterations = 1e6);

    %Variable to store the optimization results
    opt_results.best_ev = zeros(1,kmax+2);
    opt_results.best_constr = zeros(1,kmax+2);

    % Store the initial guess
    opt_results.best_ev(1) = objFunction(initGuess);
    opt_results.best_constr(1) = constraintFunction(initGuess);
    
    %First local optimization
    [opt_solution, fout] = fmincon(objFunction,initGuess,[],[],[],[],infBounds,supBounds,constraintFunction,opt_options);
    
    %Monotonic basin hopping loop
    j = 0;
    rep = 0;
    current_best = opt_solution;
    current_best_value = fout;
    current_best_constraint = constraintFunction(current_best);
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
        [opt_solution, fout] = fmincon(objFunction,new_guess,[],[],[],[],infBounds,supBounds,constraintFunction,opt_options);
    
        %Update counter of local runs
        j = j + 1;
    
        %If a better solution is found, update the result
        if current_best_constraint<=1e-3
            if fout < current_best_value && constraintFunction(opt_solution)<=1e-3
                current_best = opt_solution;
                current_best_value = fout;
                current_best_constraint = constraintFunction(opt_solution);
                rep = 0;
            else
                rep = rep + 1;
            end
        else
            if constraintFunction(opt_solution) < current_best_constraint
                current_best = opt_solution;
                current_best_value = fout;
                current_best_constraint = constraintFunction(opt_solution);
                rep = 0;
            else
                rep = rep + 1;
            end
        end
    
        %Save current best
        opt_results.best_ev(j+2) = current_best_value;
        opt_results.best_constr(j+2) = current_best_constraint;
    
        %After the specified runs, if no improvement is found, stop the algorithm
        if rep == max_rep
            break
        end
    end
    
    %Store the results of the optimization
    opt_results.local_runs = j+1;
    opt_results.repeats = rep;

end