function newsolution = removeAsteroids(solution, indexVector)

    newsolution = solution;
    
    % Validate index input
    indexVector = unique(indexVector);  % Eliminate duplicates
    indexVector(indexVector > numel(solution.Asteroid_Seq_ID)) = [];  % Remove out of range index
    
    % Remove fields
    newsolution.Asteroid_Seq_ID(indexVector) = [];
    newsolution.Flyby_Dates(indexVector) = [];
    newsolution.SC_Mass(indexVector) = [];
end