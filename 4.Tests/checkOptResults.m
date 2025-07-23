counter = 0;
arcNames = fieldnames(opt_results);
narcs = length(arcNames);

for i = 1:narcs

    if opt_results.(arcNames{i}).gaActivated == 1

        counter = counter + 1;

    end

end