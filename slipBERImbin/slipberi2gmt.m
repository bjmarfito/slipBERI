%% Code to convert slipBERI models into GMT format
% Usage: slipberi2gmt(filename) 
% Author: Bryan Marfito, 30 September 2023

function[] = slipberi2gmt(fileName)
    % Load only the necessary variables from slipBERI output file
    load(fileName, "d_InSAR", "G", "patch_mean", "faults", "locs_InSAR_latlong", "utmepicenter", "rake_mean")

    % Clear unnecessary variables
    %clearvars -except d_InSAR G patch_mean faults locs_InSAR_latlong utmepicenter rake_mean

    % Extract the x and y coordinates of the fault
    xCoord = faults(1,:);
    yCoord = faults(2,:);

    % Calculate the distance along the strike of the fault
    crossSectionAlongStrikeCoord = (xCoord.^2 + yCoord.^2) .^0.5;

    %Epicenter plotting for depth
    utmepicenterToDepth = (utmepicenter(1,2).^2 + utmepicenter(1,1).^2) .^0.5;
    utmepicenterToDepth = (utmepicenterToDepth - crossSectionAlongStrikeCoord(1));
    utmepicenterToDepth = -1.* utmepicenterToDepth./1000;
    utmepicenterToDepth = [utmepicenterToDepth -1.*utmepicenter(1,3)];
    writematrix(utmepicenterToDepth, 'hypocenter.txt', 'Delimiter', 'space');

    % Reference the along-strike distance to the minimum distance
    crossSectionAlongStrikeCoord = crossSectionAlongStrikeCoord - crossSectionAlongStrikeCoord(1);
    crossSectionAlongStrikeCoord(crossSectionAlongStrikeCoord ~= 0) = -crossSectionAlongStrikeCoord(crossSectionAlongStrikeCoord ~= 0);
    crossSectionAlongStrikeCoord = crossSectionAlongStrikeCoord';
    

    % Convert the distance along strike to km
    crossSectionAlongStrikeCoord = crossSectionAlongStrikeCoord ./1000;

    % Creates a file referenced to the hypocenter location
    % Disabled for now
    %crossSectionAlongStrikeCoordHypocenterReferenced = crossSectionAlongStrikeCoord + utmepicenterToDepth(1,1);
    
    % Extract the along-dip value and slip keep/save value for each point
    % Convert along-dip to km
    alongDipCoord = faults(8,:)' ./1000;
    alongDipCoord(alongDipCoord ~=0) = -alongDipCoord(alongDipCoord ~=0);

    % Convert slip to cm
    slip_keep_save = faults(6,:)' .*100;

    % Combine the values into a matrix for output
    gmtValueFormat = [crossSectionAlongStrikeCoord alongDipCoord slip_keep_save];

    % Write the output to a file
    outputFileName = 'slip.txt';
    writematrix(gmtValueFormat, outputFileName, 'Delimiter', 'space');

    % Combine the hypocenter referenced values into a matrix for output
    % Disabled for now
    %crossSectionAlongStrikeCoordHypocenterReferenced = -1.*crossSectionAlongStrikeCoordHypocenterReferenced;
    %gmtValueFormatHypocenterReferenced = [crossSectionAlongStrikeCoordHypocenterReferenced alongDipCoord slip_keep_save];

    % Write the output to a file
    % Disabled for now
    %outputFileNameHypocenterReferenced = 'slip_hypocenter_referenced.txt';
    %writematrix(gmtValueFormatHypocenterReferenced, outputFileNameHypocenterReferenced, 'Delimiter', 'space');

    % Slip vector
    slipVector =[crossSectionAlongStrikeCoord alongDipCoord rake_mean slip_keep_save];
    writematrix(slipVector, 'slip_vector.txt', 'Delimiter', 'space'); 

    % Generate metadata for subsurface slip
    metaDataSlip = [min(crossSectionAlongStrikeCoord) max(crossSectionAlongStrikeCoord); min(alongDipCoord) max(alongDipCoord)];
    writematrix(metaDataSlip, 'metadata_slip.txt', 'Delimiter', 'space');

    % Generate metadata for hypocenter referenced subsurface slip
    % Disabled for now
    %metaDataSlip = [min(crossSectionAlongStrikeCoordHypocenterReferenced) max(crossSectionAlongStrikeCoordHypocenterReferenced); min(alongDipCoord) max(alongDipCoord)];
    %writematrix(metaDataSlip, 'metadata_slip_hypocenter_referenced.txt', 'Delimiter', 'space');

    % Generate data, model and residual files for GMT
    lonlat = locs_InSAR_latlong(1:2,:)';
    insarData = [lonlat  d_InSAR'];

    % Surface displacement model, dhat  = Gs
    dhat = G * patch_mean;
    insarModel = [lonlat  dhat];

    % Residuals
    res = d_InSAR' - dhat;
    insarResidual = [lonlat  res];

    % Data, model and residual files for GMT in cm
    insarDataCm = [insarData(:,1:2) insarData(:,3) * 100];
    insarModelCm = [insarModel(:,1:2) insarModel(:,3) * 100];
    insarResidualCm = [insarModel(:,1:2) insarResidual(:,3) * 100];


    % Write the output to data, model, and residual files
    outputFileName = 'insar';
    writematrix(insarData, [outputFileName '_data.txt'], 'Delimiter', 'space');
    writematrix(insarModel, [outputFileName '_model.txt'], 'Delimiter', 'space');
    writematrix(insarResidual, [outputFileName '_residual.txt'], 'Delimiter', 'space');
    writematrix(insarDataCm, [outputFileName '_data_cm.txt'], 'Delimiter', 'space');
    writematrix(insarModelCm, [outputFileName '_model_cm.txt'], 'Delimiter', 'space');
    writematrix(insarResidualCm, [outputFileName '_residual_cm.txt'], 'Delimiter', 'space');

end