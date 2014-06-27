%--------------------------------------------------------------------------
% Function    : getPhiAt
% Description : This function calculates the value of the phi (father
%               wavelet) at given points by using the gridded interpolant.
% 
% Inputs  :
%     - points  : points at which the value of phi is to be calculated.
% 
% Outputs :
%     - phiValue : the value of phi at the given points.
% 
% Usage   : This function is used to calculate the value of phi at given
%           points using the gridded interpolant.
%           
%           This function uses the static (global) DensityVars.
%          
% Author(s):
%   - Gedeon Nyengele (gnyengele3@mail.gatech.edu)
%
% Date Published: 6/26/2014
%-------------------------------------------------------------------------- 
function phiValue = getPhiAt( points )
    
% Access the static (global) DensityVars.
global DensityVars;

% Calculate phi values by interpolation.
phiValue = DensityVars.PhiInterpolant( points );

% Get range of allowable domain values.
range = DensityVars.WaveSupport;

% Correct false interpolated values of phi.
phiValue( points < range(1) ) = 0;
phiValue( points > range(2) ) = 0;

end % end function getPhiAt