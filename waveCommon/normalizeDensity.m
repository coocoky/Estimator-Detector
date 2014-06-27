%--------------------------------------------------------------------------
% Function:    normalizeDensity
% Description: Ensures the density integrates to 1.
%
% Inputs:
%   densDom     - Vector with the limit points of the density domain.
%   densEst     - This refers to the density values that are going to be 
%                    normalized.
%   dt             - This is the discretization value in this case.    
% Outputs:
%   normDens       - Vector of normalized density function values.
% Usage: This function normalizes the density function values to ensure
%        that they are non-negative and integrate to 1.
%
% Authors(s):
%   Eddy Ihou(ihouk2002@my.fit.edu)
%   Mark Moyou(mmoyou@my.fit.edu)
%--------------------------------------------------------------------------
function normDens = normalizeDensity(densDom,densEst,dt)

% Set first iteration value equal to the initial estimate and k = 0.
normDens = densEst;
iter = 0;
threshold = 10e-16;
% Density domain spread.
domainSize = densDom(2) - densDom(1);

% Loop until the maximum number of iterations reached
while (iter < 1000)
    iter = iter + 1;
    
    % Set all negative values to zero
    normDens = max(normDens,0);
    % Check to see if density integrates to 1
    densityIntegral = abs(sum(normDens)*dt - 1);

    if (densityIntegral < threshold)
        % If density estimate == 1 then stop and return density.  
        return;
    else 
        intgDens = sum(normDens)*dt;  
        % Else set new estimate to previous manipulated estimate - (Ck+1 -
        % 1)/weighting function.
        normDens = normDens - ((intgDens - 1) / domainSize); 
    end
end

end