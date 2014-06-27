%--------------------------------------------------------------------------
% Function:    calculateDensity
% Description: Calculates the density with the updated coefficients.
%
% Inputs:
%   densityDomain     - Maximum and minimum range of the samples. Samples
%                       outside this range will be discarded.
%   wavelet           - Name of wavelet to use for density approximation.
%                       Use matlab naming convention for wavelets.
%   startLevel        - Starting level for the the father wavelet
%                       (i.e. scaling function).  
%   stopLevel         - Last level for mother wavelet scaling. The start
%                       level is same as the father wavelet's.
%   coefficients      - Vector of coefficients.  
%   waveletFlag       - If this flag is on then density estimation is done     
%                       with scaling + wavelets. The default is density
%                       estimation with the scaling function only. 
%                       waveletFlag = 1; Wavelet is on. 
%                       waveletFlag = 0; Wavelet is off.
%   discretization    - Discretization of the domain.
%   scalTranslates    - Scaling translates.
%   waveTranslates    - Wavelet translates.
%   subSample         - This is a vector of samples that the user inputs in 
%                       order to obtain p(x) values.   
% Outputs:
%   density           - Values of the density function on domain.
%   subDensity        - The p(x) values that correspond to the subSample.
% Usage:
%
% Authors(s):
%   Eddy Ihou(ihouk2002@my.fit.edu)
%   Mark Moyou(mmoyou@my.fit.edu)
%   Gedeon Nyengele (gnyengele3@mail.gatech.edu)
%--------------------------------------------------------------------------
function density = calculateDensity()

global DensityVars;
persistent savedPhis;
persistent savedPsis;

% Create a vector for the density of the domain.
domain = (DensityVars.DensityDomain(1) : DensityVars.Discretization : DensityVars.DensityDomain(2));

% Set up the probability density function as 0 over the domain.
density = zeros(1, length(domain));
lengthDomain = length(domain);


    
% SAVE PHIs
% This is ran only once during the entire density estimation process.
if( isempty(savedPhis) )
    pts = repmat( domain',1, length(DensityVars.ScalTranslates) );
    kVals = repmat( DensityVars.ScalTranslates, lengthDomain, 1 );
    xVals = 2.^DensityVars.StartLevel*pts - kVals;
    phis = getPhiAt( xVals(:) );
    savedPhis = 2.^(DensityVars.StartLevel / 2) * reshape(phis, size(xVals));        

end % if( isempty(savedPhis) )  ----- END SAVE PHIs.

% Calculate density using father wavelet.
coeffs = repmat(DensityVars.Coefficients{1,1}, lengthDomain, 1);
dens = sum( coeffs .* savedPhis, 2)';

if( DensityVars.WaveletFlag == 1 ) % Wavelet is ON
    
    % SAVE PSIs.
    % This is ran only once during the entire density estimation process.
    if( isempty(savedPsis) )
        numLevels = 1 + DensityVars.StopLevel - DensityVars.StartLevel;
        savedPsis = cell(numLevels, 1);
        
        for J = DensityVars.StartLevel : DensityVars.StopLevel
            pos = 1 + J - DensityVars.StartLevel;
            pts = repmat( domain',1, length(DensityVars.WaveTranslates{pos,1}) );
            kVals = repmat( DensityVars.WaveTranslates{pos,1}, lengthDomain, 1 );
            xVals = 2.^J*pts - kVals;
            psis = getPsiAt( xVals(:) );
            savedPsis{pos} = 2.^(J / 2) * reshape(psis, size(xVals));
        end % for l = DensityVars.StartLevel : DensityVars.StopLevel
        
        
    end % if( isempty(savedPsis) ) ----- END SAVE PSIs.
    
    % Compute addendum density.
    for J = DensityVars.StartLevel : DensityVars.StopLevel
        pos = 1 + J - DensityVars.StartLevel;
        coeffs = repmat(DensityVars.Coefficients{pos + 1,1}, lengthDomain, 1);
        phis = savedPsis{pos, 1};
        subdens = sum( coeffs .* phis, 2 )';
        dens = dens + subdens;
    end
    
end % if( DensityVars.WaveletFlag == 1 )

% Normalize the density
density = normalizeDensity( DensityVars.DensityDomain,...
                            dens,...
                            DensityVars.Discretization );


end % function density = calculateDensity()