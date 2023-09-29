%% Code to implement burn-in for the results from slipBERI
% Usage: burnin_results(filename, no. of burn-in)
% Author: Bryan Marfito, Sept 29, 2023

function [] = burnin_results(fileName, noOfsamplesToBurnIn)
    % Load the data from the specified file
    try
        load(fileName)
    catch exception
        error('Error loading file: %s', exception.message)
    end
    
    % Set the number of burn-in samples to remove
    burnInRemoveNumber = noOfsamplesToBurnIn;
    
    % Print a message indicating how many burn-in samples will be removed
    fprintf('Removing %d burn-in samples from %s\n', burnInRemoveNumber, fileName)
    
    % Call the remove_burn_in function to remove the specified number of burn-in samples
    try
        remove_burn_in
    catch exception
        error('Error removing burn-in samples: %s', exception.message)
    end
end