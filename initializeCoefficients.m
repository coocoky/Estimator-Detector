
function  initializeCoefficients()

global DensityVars;

% Set all coefficients to zero.                                                        
DensityVars.Coefficients = cell((DensityVars.StopLevel - DensityVars.StartLevel + 2),2); % Creating coefficients cell.                                                        
DensityVars.Coefficients{1,1} = zeros(size(DensityVars.ScalTranslates)); % Setting scaling coefficients to zero.
DensityVars.Coefficients{1,2} = DensityVars.StartLevel;

if (DensityVars.WaveletFlag == 1) % Scaling and wavelets.
    for j = DensityVars.StartLevel : DensityVars.StopLevel        
        coeffAtResj = DensityVars.WaveTranslates{((j - DensityVars.StartLevel) + 1),1};
        coeffIndVal = ((j - DensityVars.StartLevel) + 2) ; % Only wavelet coefficients
        DensityVars.Coefficients{coeffIndVal,1} = zeros(1,length(coeffAtResj));
        DensityVars.Coefficients{coeffIndVal,2} = j;
    end
end

end % initializingTheCoefficients()