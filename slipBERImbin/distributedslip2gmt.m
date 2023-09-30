%% Code to convert slipBERI slip models into GMT format
% Usage: slip2gmt(filename) converts a slipberi file to a GMT subsurface file
% Author: Bryan Marfito, 30 September 2023

function[] = distributedslip2gmt(fileName)
    % Load the slipBERI file
    load(fileName)

    % Clear unnecessary variables
    clearvars -except d_InSAR G patch_mean faults locs_InSAR_latlong

    % Extract the x and y coordinates of the fault
    xCoord = faults(1,:);
    yCoord = faults(2,:);

    % Calculate the distance along the strike of the fault
    crossSectionAlongStrikeCoord = (xCoord.^2 + yCoord.^2) .^0.5;
    crossSectionAlongStrikeCoord = crossSectionAlongStrikeCoord - crossSectionAlongStrikeCoord(1);
    crossSectionAlongStrikeCoord(crossSectionAlongStrikeCoord ~= 0) = -crossSectionAlongStrikeCoord(crossSectionAlongStrikeCoord ~= 0);
    crossSectionAlongStrikeCoord = crossSectionAlongStrikeCoord';

    % Extract the along-dip value and slip keep/save value for each point
    alongDipCoord = faults(8,:)';
    slip_keep_save = faults(6,:)';

    % Combine the values into a matrix for output
    gmtValueFormat = [crossSectionAlongStrikeCoord alongDipCoord slip_keep_save];

    % Write the output to a file
    writematrix(gmtValueFormat, 'slip.txt', 'Delimiter', 'space');

    % Generate metadata for subsurface slip
    meteDataSlip = [min(crossSectionAlongStrikeCoord) max(crossSectionAlongStrikeCoord); min(alongDipCoord) max(alongDipCoord)];
    writematrix(meteDataSlip, 'metadata_slip.txt', 'Delimiter', 'space');

    % Extract only lon and lat values
    lonlat = transpose(locs_InSAR_latlong(1:2,:));
    
    % Generate data, model and residual files for GMT

    insarData = [lonlat  d_InSAR'];

    % Surface displacement model, dhat  = Gs
    dhat = G * patch_mean;
    insarModel = [lonlat  dhat];

    % Residuals
    res = d_InSAR' - dhat;
    insarResidual = [lonlat  res];

    % Data, model and residual files for GMT in cm
    insarDataCm = [insarData(:,1) insarData(:,2 ) insarData(:,3) * 100];
    insarModelCm = [insarModel(:,1) insarModel(:,2) insarModel(:,3) * 100];
    insarResidualCm = [insarModel(:,1) insarModel(:,2) insarResidual(:,3) * 100];

    % Write the output to data, model, and residual files
    writematrix(insarData, 'insar_data.txt', 'Delimiter', 'space');
    writematrix(insarModel, 'insar_model.txt', 'Delimiter', 'space');
    writematrix(insarResidual, 'insar_residual.txt', 'Delimiter', 'space');
    writematrix(insarDataCm, 'insar_data_cm.txt', 'Delimiter', 'space');
    writematrix(insarModelCm, 'insar_model_cm.txt', 'Delimiter', 'space');
    writematrix(insarResidualCm, 'insar_residual_cm.txt', 'Delimiter', 'space');

end