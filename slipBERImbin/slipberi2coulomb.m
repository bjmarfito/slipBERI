function[] = slipberi2coulomb(fileName)
    load(fileName,"disloc_model","rake_mean","patch_mean")

    % Retain disloc_model, rake_mean, and slip_mean
    faults = disloc_model';
    noSubPatches = size(faults);

    for faultPatch = 1:noSubPatches

        xCenter = faults(faultPatch,1);
        yCenter = faults(faultPatch,2);
        faultLength = faults(faultPatch,7);
        faultStrike = deg2rad(faults(faultPatch,3));
        faultDip = deg2rad(faults(faultPatch,4));
        faultTop = faults(faultPatch,8);
        faultBottom = faults(faultPatch,9);

        xPoint1 = xCenter - 0.5*faultLength*sin(faultStrike);
        xPoint2 = xCenter + 0.5*faultLength*sin(faultStrike);
        yPoint1 = yCenter - 0.5*faultLength*cos(faultStrike);
        yPoint2 = yCenter + 0.5*faultLength*cos(faultStrike);

        xTop1 = xPoint1 + faultTop*cos(faultStrike)/tan(faultDip);
        xTop2 = xPoint2 + faultTop*cos(faultStrike)/tan(faultDip);
        yTop1 = yPoint1 - faultTop*sin(faultStrike)/tan(faultDip);
        yTop2 = yPoint2 - faultTop*sin(faultStrike)/tan(faultDip);

        xBot1 = xPoint1 + faultBottom *cos(faultStrike)/tan(faultDip);
        xBot2 = xPoint2 + faultBottom*cos(faultStrike)/tan(faultDip);
        yBot1 = yPoint1 - faultBottom*sin(faultStrike)/tan(faultDip);
        yBot2 = yPoint2 - faultBottom*sin(faultStrike)/tan(faultDip);

        patchX(1:4,faultPatch) = [xTop1 xTop2 xBot2 xBot1]'/1000;
        patchY(1:4,faultPatch) = [yTop1 yTop2 yBot2 yBot1]'/1000;
        patchZ(1:4,faultPatch) = [-faultTop  -faultTop  -faultBottom  -faultBottom]'/1000;

    end
end