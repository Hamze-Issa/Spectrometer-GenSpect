function countsVector = radianceToCounts(path, pixelWidth, solidAngle, lensArea, intTime)

    convolvedRadiance = instconv(path.radiance, l2nu(1300) - l2nu(900), l2nu(1300) - l2nu(900));
    countsVector = convolvedRadiance .* (pixelWidth * solidAngle * lensArea * intTime);
    
end