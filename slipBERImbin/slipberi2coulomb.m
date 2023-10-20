function[] = slipberi2coulomb(fileName)
    % Load only the necessary variables from slipBERI output file
    load(fileName,"disloc_model","rake_mean","patch_mean")

    faults = disloc_model';
    [noSubPatches,~] = size(faults);

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

        % Calculate the coordinates of the top of the fault patch using the Aki-Richards coordinates
        xTop1 = xPoint1 + faultTop*cos(faultStrike)/tan(faultDip);
        xTop2 = xPoint2 + faultTop*cos(faultStrike)/tan(faultDip);
        yTop1 = yPoint1 - faultTop*sin(faultStrike)/tan(faultDip);
        yTop2 = yPoint2 - faultTop*sin(faultStrike)/tan(faultDip);

    
        % Calculate the coordinates of the bottom of the fault patch using the Aki-Richards coordinates
        xBot1 = xPoint1 + faultBottom*cos(faultStrike)/tan(faultDip);
        xBot2 = xPoint2 + faultBottom*cos(faultStrike)/tan(faultDip);
        yBot1 = yPoint1 - faultBottom*sin(faultStrike)/tan(faultDip);
        yBot2 = yPoint2 - faultBottom*sin(faultStrike)/tan(faultDip);

        % Stores the coordinates of the fault per x-, y-, and z- coordinates
        % and converts them from meters to kilometers
        patchX(1:4,faultPatch) = [xTop1 xTop2 xBot2 xBot1]'/1000;
        patchY(1:4,faultPatch) = [yTop1 yTop2 yBot2 yBot1]'/1000;
        patchZ(1:4,faultPatch) = [-faultTop  -faultTop  -faultBottom  -faultBottom]'/1000;

        % Kode value for subfaults
        kode = ones(noSubPatches,1).*100;

        % Dip value for subfaults
        faultDipPatches = ones(noSubPatches,1).*rad2deg(faultDip);

    end

    % Save the fault patches in a .inr file
    matrix2inr = [patchX(:,1)' patchY(:,1)' patchX(:,2)' patchY(:,2)' kode(1,:) rake_mean patch_mean faultDipPatches patchZ(:,1)' patchZ(:,3)'];
    save -ascii slipBERI2Coulomb.inr matrix2inr
end