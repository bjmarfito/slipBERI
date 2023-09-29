%% Code to convert slipBERI slip models into GMT format
% Usage: slip2gmt(filename) converts a slipberi file to a GMT subsurface file

function[] = distributedslip2gmt(fileName)
    % Load the slipBERI file
    load(fileName)

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

    % Generate metadat for subsurface slip
    meteDataSlip = [min(crossSectionAlongStrikeCoord) max(crossSectionAlongStrikeCoord); min(alongDipCoord) max(alongDipCoord)];
    writematrix(meteDataSlip, 'metadata_slip.txt', 'Delimiter', 'space');
end