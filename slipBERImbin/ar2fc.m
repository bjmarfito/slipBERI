% Function to transform Aki-Richards to fault centered coordinates
function[lonc, latc] = ar2fc(lona2r, lata2r, faultStrike, faultLength)

    arcLen = km2deg(0.5*faultLength);
    rotateCenterToAR = faultStrike + 180;
    [lonc, latc] = reckon(lata2r, lona2r, arcLen, rotateCenterToAR);
    
end

    
