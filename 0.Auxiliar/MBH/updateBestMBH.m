function [newBest,newBest_value,newBest_constraint] = updateBestMBH(newSolution,current_best,objFunction,constraintFunction)

    % Get the value and constraint of new solution
    newSolution_value = objFunction(newSolution);
    newSolution_constraint = max(constraintFunction(newSolution));

    % Get the value and constraint of current best solution
    current_best_value = objFunction(current_best);
    current_best_constraint = max(constraintFunction(current_best));

    % Apply the update logic
    if current_best_constraint<=1e-3 && newSolution_value < current_best_value && newSolution_constraint<=1e-3
        newBest = newSolution;
        newBest_value = newSolution_value;
        newBest_constraint = newSolution_constraint;
    elseif current_best_constraint>1e-3 && newSolution_constraint < current_best_constraint
        newBest = newSolution;
        newBest_value = newSolution_value;
        newBest_constraint = newSolution_constraint;
    else
        newBest = current_best;
        newBest_value = current_best_value;
        newBest_constraint = current_best_constraint;
    end
end