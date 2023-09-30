%% Code to implement burn-in for the results from slipBERI
% Usage: burnin_results(filename, no. of burn-in)
% Author: Bryan Marfito, 29 Sept 2023

function [] = burnin_results(fileName, noOfsamplesToBurnIn)
    % Load the data from the specified file
    try
        load(fileName)
    catch exception
        error('Error loading file: %s', exception.message)
    end

    % Close any opened figures to save memory
    close all
    
    % Set the number of burn-in samples to remove
    burn_in_remove_number = noOfsamplesToBurnIn;
    
    % Print a message indicating how many burn-in samples will be removed
    fprintf('Removing %d burn-in samples from %s .\n', burn_in_remove_number, fileName)
    
    % Call the remove_burn_in function to remove the specified number of burn-in samples
    try
        remove_burn_in
    catch exception
        error('Error removing burn-in samples: %s', exception.message)
    end

    % Save slip_keep results into a text file
    save slip_keep.txt slip_keep -ascii

end