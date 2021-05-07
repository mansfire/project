function newDat = spectrumInterpolation(data, Fs, Fl, neighborsToSample, neighborsToReplace)
    spectrum = fft(data);
    mag = abs(spectrum); % We use the real spectrum to interpolate and remove the powerline noise.
    phase = angle(spectrum);
    
    binStepSize = length(mag) / Fs;
    nextPowerBin = binStepSize * Fl;
    neighborsToSample = uint32(binStepSize * neighborsToSample);
    neighborsToReplace =  uint32(binStepSize * neighborsToReplace);
    % We have to take care of each end of the spectrum differently.
    % For the first half, the spectra is not reversed.    
    nearestPowerlineHarmonic = nextPowerBin;
    binStart = nearestPowerlineHarmonic;
    
    for i=binStart:nextPowerBin:length(mag)/2
        neighborSamples = mag(i-neighborsToSample:i+neighborsToSample, :);
        neighborhoodAverage = median(neighborSamples);
        for j = 1:width(data)
            mag(i-neighborsToReplace:i+neighborsToReplace, j) = neighborhoodAverage(j);
        end
    end
    
    nyquistFrequency = Fs / 2;
    nearestPowerlineHarmonic = mod(nyquistFrequency, Fl); % This gives us the distance to the nyquist frequency starting from the middle.
    binStart = (length(mag)/2) + (nearestPowerlineHarmonic * binStepSize);
    
    % Now replace the other side

    for i=binStart:nextPowerBin:(length(mag) - nextPowerBin)
        neighborSamples = mag(i-neighborsToSample:i+neighborsToSample, :);
        neighborhoodAverage = median(neighborSamples);
        for j = 1:width(data)
            mag(i-neighborsToReplace:i+neighborsToReplace, j) = neighborhoodAverage(j);
        end
    end

    
    % Create a new signal based on euler's formula.    
    newDat = mag.*exp(1i.*phase);
    newDat = real(ifft(newDat));
end