% Function to transform Aki-Richards to fault centered coordinates
function[latc, lonc] = ar2fc(lata2r, lona2r, faultStrike, faultLength)

    arcLen = km2deg(0.5*faultLength);
    [latc, lonc] = reckon(lata2r, lona2r, arcLen, faultStrike);

end

    
