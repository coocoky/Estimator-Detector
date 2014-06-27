%--------------------------------------------------------------------------
% Function    : getPsiAt
% Description : This function calculates the value of the psi (mother
%               wavelet) at given points by using the gridded interpolant.
% 
% Inputs  :
%     - points  : points at which the value of psi is to be calculated.
% 
% Outputs :
%     - psiValue : the value of psi at the given points.
% 
% Usage   : This function is used to calculate the value of psi at given
%           points using the gridded interpolant.
%           
%           This function uses the static (global) DensityVars.
%          
% Author(s):
%   - Gedeon Nyengele (gnyengele3@mail.gatech.edu)
%
% Date Published: 6/26/2014
%-------------------------------------------------------------------------- 
function psiValue = getPsiAt( points )
    
% Access the static (global) DensityVars.
global DensityVars;

% Calculate phi values by interpolation.
psiValue = DensityVars.PsiInterpolant( points );

% Get range of allowable domain values.
range = DensityVars.WaveSupport;

% Correct false interpolated values of phi.
psiValue( points < range(1) ) = 0;
psiValue( points > range(2) ) = 0;

end % end function getPhiAt