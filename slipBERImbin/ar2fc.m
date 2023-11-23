% Function to transform Aki-Richards to fault centered coordinates
function[latc, lonc] = ar2fc(lata2r, lona2r, faultStrike, faultLength)

    arcLen = km2deg(0.5*faultLength);
    rotateCenterToAR = faultStrike + 180;
    [latc, lonc] = reckon(lata2r, lona2r, arcLen, rotateCenterToAR);
    
end

    
