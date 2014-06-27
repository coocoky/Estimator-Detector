function getTranslates()

global DensityVars;

% Find translation range dependent on the domain of the density and the
% wavelet basis.
transRange = translationRange(DensityVars.StartLevel);
DensityVars.ScalTranslates = transRange(1) : transRange(2);
DensityVars.WaveTranslates = [];

if (DensityVars.WaveletFlag == 1) % wavelet is on. 
    % Creating a cell to hold wavelet translates. 
    numWavelets = (DensityVars.StopLevel - DensityVars.StartLevel) + 1;
    DensityVars.WaveTranslates = cell(numWavelets,2);
    for i = DensityVars.StartLevel : DensityVars.StopLevel
        waveTransRange = translationRange(i);
        % Storing the wavelet translate values. 
        DensityVars.WaveTranslates {((i - DensityVars.StartLevel) + 1),1} = waveTransRange(1) : waveTransRange(2);
        % Storing the corresponding resolution levels. 
        DensityVars.WaveTranslates {((i - DensityVars.StartLevel) + 1),2} = i;
    end % i = startLevel : stopLevel

end % (waveletFlag == 1) % wavelet is on.  


end % function [scalTranslates waveTranslates] = getTranslates(densityDomain,...