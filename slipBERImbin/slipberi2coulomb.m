%% Code to convert slipBERI output file to Coulomb input file
% slipBERI output file is in the form of a .mat file and generates a Coulomb input file in .inr format
% Usage: slipberi2coulomb(name of slipBERI output file)
% Example: slipberi2coulomb('slipBERI_output.mat')
% Author: Bryan Marfito, 11 Nov 2023

function[] = slipberi2coulomb(fileName)
    % Load only the necessary variables from the slipBERI output file
    load(fileName,"disloc_model","rake_mean","patch_mean", "utmepicenter", "EQ_epicenter")

    % Incremet of the grid in degrees
    distanceIncrement = 0.05;

    %Create a +- 1 degree extent from the epicenter for the Coulomb stress transfer models
    minGrid = EQ_epicenter(1:2) - 1;
    maxGrid = EQ_epicenter(1:2) + 1;

    % Convert minimum and maximum extent from geocoordinates to UTM coordinates
    [minGridLonUTM, minGridLatUTM, ~] = ll2utm(minGrid(2), minGrid(1));
    [maxGridLonUTM, maxGridLatUTM, ~] = ll2utm(maxGrid(2), maxGrid(1));

    % Calculate the number of grid points in the Caretsian coordinates (x- and y-directions)
    % and uses the epicenter as the reference point
    lonIncrement  =  EQ_epicenter(1) + distanceIncrement;
    latIncrement  =  EQ_epicenter(2) + distanceIncrement;
    [lonIncrementUTM, latIncrementUTM, ~] = ll2utm(latIncrement, lonIncrement);
    lonIncrementUTM = lonIncrementUTM - utmepicenter(1);
    latIncrementUTM = latIncrementUTM - utmepicenter(2);
    minGridLonUTM = minGridLonUTM - utmepicenter(1);
    minGridLatUTM = minGridLatUTM - utmepicenter(2);
    maxGridLonUTM = maxGridLonUTM - utmepicenter(1);
    maxGridLatUTM = maxGridLatUTM - utmepicenter(2);
    
    % Array for grid parameters in kilometers which will be needed by Coulomb for stress transfer calculation
    gridArray = [minGridLonUTM; minGridLatUTM; maxGridLonUTM; maxGridLatUTM; lonIncrementUTM; latIncrementUTM] ./1000;

    % Poisson's ratio
    nu = 0.25;

    % Plotting parameters
    plotParameters = [2; 1; 10000];

    % Young's modulus in Pa
    eta = 800000;

    % Depth in km where the stress is calculated
    % This parameter can be changed later on in Coulomb software
    depthKmStressCalc = 4.0;

    % Friction coefficient
    frictionCoefficient = 0.6;

    % Regional stress field(3x4) matrix, default values from Coulomb 3.3 software
    % Modified from global_variable_explanation.m 
    % in the Coulomb 3.3 software
    % Useful only when stress transfer is calculated in optimally-oriented faults
    % The values can be changed later on in Coulomb software
    % regionalStress(1,1): orientation of sigma-1 (degree)
    % regionalStress(1,2): plunge of sigma-1 (degree)
    % regionalStress(1,3): magnitude at the surface (bar)
    % regionalStress(1,4): depth gradient (S0 + gradient * depth (km))
    % regionalStress(2,1): orientation of sigma-2 (degree)
    % regionalStress(2,2): plunge of sigma-2 (degree)
    % regionalStress(2,3): magnitude at the surface (bar)
    % regionalStress(2,4): depth gradient (S0 + gradient * depth (km))
    % regionalStress(3,1): orientation of sigma-3 (degree)
    % regionalStress(3,2): plunge of sigma-3 (degree)
    % regionalStress(3,3): magnitude at the surface (bar)
    % regionalStress(3,4): depth gradient (S0 + gradient * depth (km))
    regionalStress = [19 -0.01 100 0; 89.99 89.99 30 0; 109 -0.01 0 0];

    % Extract the parameters of subpatches from the slipBERI output file
    faults = disloc_model';
    [noSubPatches,~] = size(faults);

    % Initialize the patch coordinates to preallocate the no of arrays needed for
    %storing 
    patchX = zeros(4,noSubPatches);
    patchY = zeros(4,noSubPatches);
    patchZ = zeros(4,noSubPatches);

    % Create ID for the main fault to be analyzed
    faultSegmentID = ones(noSubPatches,1);

    % Generetes a default cross-section parameters
    % This can be changed in the Coulomb software
    crossSectionParams = [-16; -16; 18; 26; 1; 30; 1];

    % Loop through each fault patch and calculate the four corner coordinates of the fault subpatches
    % using the given center coordinates, strike, dip, top, and bottom of the fault.
    for faultPatch = 1:noSubPatches

        % Extract fault parameters from slipBERI output file
        % x-center, y-center, strike, dip, rake, slip, length, top, bottom
        xCenter = faults(faultPatch,1);
        yCenter = faults(faultPatch,2);
        faultLength = faults(faultPatch,7);
        faultStrike = deg2rad(faults(faultPatch,3));
        faultDip = deg2rad(faults(faultPatch,4));
        faultTop = faults(faultPatch,8);
        faultBottom = faults(faultPatch,9);

        % Calculate the Aki-Richards x- and y- coordinates of the fault patches using the fault center
        xPoint1 = xCenter - 0.5*faultLength*sin(faultStrike);
        xPoint2 = xCenter + 0.5*faultLength*sin(faultStrike);
        yPoint1 = yCenter - 0.5*faultLength*cos(faultStrike);
        yPoint2 = yCenter + 0.5*faultLength*cos(faultStrike);

        % Calculate the x- and y- coordinates of the top of each fault subpatch
        xTop1 = xPoint1 + faultTop*cos(faultStrike)/tan(faultDip);
        xTop2 = xPoint2 + faultTop*cos(faultStrike)/tan(faultDip);
        yTop1 = yPoint1 - faultTop*sin(faultStrike)/tan(faultDip);
        yTop2 = yPoint2 - faultTop*sin(faultStrike)/tan(faultDip);

        % Calculate the x- and y- coordinates coordinates of the bottom of each fault subpatch
        xBot1 = xPoint1 + faultBottom*cos(faultStrike)/tan(faultDip);
        xBot2 = xPoint2 + faultBottom*cos(faultStrike)/tan(faultDip);
        yBot1 = yPoint1 - faultBottom*sin(faultStrike)/tan(faultDip);
        yBot2 = yPoint2 - faultBottom*sin(faultStrike)/tan(faultDip);

        % Stores the coordinates of each fault subpatch per x-, y-, and z- coordinates    
        patchX(1:4,faultPatch) = [xTop1 xTop2 xBot2 xBot1]';
        patchY(1:4,faultPatch) = [yTop1 yTop2 yBot2 yBot1]';
        patchZ(1:4,faultPatch) = [-faultTop  -faultTop  -faultBottom  -faultBottom]';

        % Converts dip values for subpatches to degrees again since Coulomb uses angles in degrees
        faultDipPatches = ones(noSubPatches,1).*rad2deg(faultDip);
    end

    % Assign 100 as the kode value for faults
    kode = ones(noSubPatches,1).*100;

    % Make the epicenter as the reference point of thex- and y- coordinates 
    % and converts them from meters to kilometers
    patchX = (patchX - utmepicenter(1)) ./1000;
    patchY = (patchY - utmepicenter(2))./1000;

    % Converts the depth of the subpatches to kilometers
    patchZ = patchZ ./1000;

    % Creates a comment name for the fault segment
    faultNo = "   Fault 1";

    % Put all the fault parameters in an array
    matrix2inr = [faultSegmentID patchX(1,:)' patchY(1,:)' patchX(2,:)' patchY(2,:)' kode rake_mean patch_mean faultDipPatches -1.*patchZ(1,:)' -1.*patchZ(3,:)'];
    
    % Generate the Coulomb input file in .inr format using the format given in the Coulomb manual
    % and using the parameters specified above
    fileID = fopen('slipBERI2Coulomb.inr','wt');
    fprintf(fileID,'header line 1 \n');
    fprintf(fileID,'header line 2 \n');
    fprintf(fileID,'#reg1=  0  #reg2=  0  #fixed= %3i  sym=  1\n',noSubPatches);
    fprintf(fileID,' PR1=%12.3f     PR2=%12.3f   DEPTH=%12.3f\n',nu,nu,depthKmStressCalc);
    fprintf(fileID,'  E1=%15.3e   E2=%15.3e\n',eta,eta);
    fprintf(fileID,'XSYM=       .000     YSYM=       .000\n');
    fprintf(fileID,'FRIC=%15.3f\n',frictionCoefficient);
    fprintf(fileID,'S1DR=%15.3f',regionalStress(1,1));
    fprintf(fileID,' S1DP=%15.3f',regionalStress(1,2));
    fprintf(fileID,' S1IN=%15.3f',regionalStress(1,3));
    fprintf(fileID,' S1GD=%15.3f\n',regionalStress(1,4));
    fprintf(fileID,'S2DR=%15.3f',regionalStress(2,1));
    fprintf(fileID,' S2DP=%15.3f',regionalStress(2,2));
    fprintf(fileID,' S2IN=%15.3f',regionalStress(2,3));
    fprintf(fileID,' S2GD=%15.3f\n',regionalStress(2,4));
    fprintf(fileID,'S3DR=%15.3f',regionalStress(3,1));
    fprintf(fileID,' S3DP=%15.3f',regionalStress(3,2));
    fprintf(fileID,' S3IN=%15.3f',regionalStress(3,3));
    fprintf(fileID,' S3GD=%15.3f\n\n',regionalStress(3,4));
    fprintf(fileID,'  #   X-start    Y-start     X-fin      Y-fin   Kode   rake     netslip   dip angle     top        bot\n');
    fprintf(fileID,'xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n');
    for faultPatch = 1:noSubPatches
        fprintf(fileID, '%3i %10.4f %10.4f %10.4f %10.4f %3i %10.4f %10.4f %10.4f %10.4f %10.4f %s\n',matrix2inr(faultPatch,:), faultNo);
    end

    fprintf(fileID,'  \n');
    fprintf(fileID,'    Grid Parameters\n');
    fprintf(fileID,'  1  ----------------------------  Start-x = %16.7f\n',gridArray(1,1));
    fprintf(fileID,'  2  ----------------------------  Start-y = %16.7f\n',gridArray(2,1));
    fprintf(fileID,'  3  --------------------------   Finish-x = %16.7f\n',gridArray(3,1));
    fprintf(fileID,'  4  --------------------------   Finish-y = %16.7f\n',gridArray(4,1));
    fprintf(fileID,'  5  ------------------------  x-increment = %16.7f\n',gridArray(5,1));
    fprintf(fileID,'  6  ------------------------  y-increment = %16.7f\n',gridArray(6,1));
    fprintf(fileID,'     Size Parameters\n');
    fprintf(fileID,'  1  --------------------------  Plot size = %16.7f\n',plotParameters(1,1));
    fprintf(fileID,'  2  --------------  Shade/Color increment = %16.7f\n',plotParameters(2,1));
    fprintf(fileID,'  3  ------  Exaggeration for disp.& dist. = %16.7f\n',plotParameters(3,1));
    fprintf(fileID,'  \n');

    fprintf(fileID,'     Cross section default\n');
    fprintf(fileID,'  1  ----------------------------  Start-x = %16.7f\n',crossSectionParams(1,1));
    fprintf(fileID,'  2  ----------------------------  Start-y = %16.7f\n',crossSectionParams(2,1));
    fprintf(fileID,'  3  --------------------------   Finish-x = %16.7f\n',crossSectionParams(3,1));
    fprintf(fileID,'  4  --------------------------   Finish-y = %16.7f\n',crossSectionParams(4,1));
    fprintf(fileID,'  5  ------------------  Distant-increment = %16.7f\n',crossSectionParams(5,1));
    fprintf(fileID,'  6  ----------------------------  Z-depth = %16.7f\n',crossSectionParams(6,1));
    fprintf(fileID,'  7  ------------------------  Z-increment = %16.7f\n',crossSectionParams(7,1));

    fprintf(fileID,'     Map info\n');
    fprintf(fileID,'  1  ---------------------------- min. lon = %16.7f\n',minGrid(1));
    fprintf(fileID,'  2  ---------------------------- max. lon = %16.7f\n',maxGrid(1));
    fprintf(fileID,'  3  ---------------------------- zero lon = %16.7f\n',EQ_epicenter(1));
    fprintf(fileID,'  4  ---------------------------- min. lat = %16.7f\n',minGrid(2));
    fprintf(fileID,'  5  ---------------------------- max. lat = %16.7f\n',maxGrid(2));
    fprintf(fileID,'  6  ---------------------------- zero lat = %16.7f\n',EQ_epicenter(2));

    fclose(fileID);
end