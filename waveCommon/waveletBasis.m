%--------------------------------------------------------------------------
% Function:    waveletBasis
% Description: Evaluates the input value using the given basis.
%
% Inputs:
%   x           - input value to be evaluated at the given basis.
%   parent      - either 'father' or 'mother' to indicate whether we want
%                 the scaling function or wavelet function value.
%   wName       - 3 to 4 charater code name of wavelet for density 
%                 approximation.
%   delta       - Discretization parameter.
% Outputs:
%   f           - output value f(x) for the given basis.
%
% Usage:
%
% Authors(s):
%   Adrian M. Peter
%
% Reference:
% A. Peter and A. Rangarajan, “Maximum likelihood wavelet density estimation 
% with applications to image and shape matching,” IEEE Trans. Image Proc., 
% vol. 17, no. 4, pp. 458–468, April 2008.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Copyright (C) 2009 Adrian M. Peter (adrian.peter@gmail.com)
%
%     This file is part of the WDE package.
%
%     The source code is provided under the terms of the GNU General 
%     Public License as published by the Free Software Foundation version 2 
%     of the License.
%
%     WDE package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with WDE package; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
%     USA
%--------------------------------------------------------------------------
function f = waveletBasis(x, parent, wName)

% Flags
persistent moth; 
persistent fath; 

% Permanently holds the table of values for the scaling and wavelet
% functions.
persistent phi;
persistent psi;
persistent supp;

% Only need to calculate the values once.
if( isempty(fath) & isempty(moth) )
    fath = 1; moth = 1;
    load([wName 'Tables']);
end %if( (fath==0) &(moth==0) )

% Calculate the basis values.
if(strcmp(lower(parent),'father'))
    % Interpolate and get function values for input x values.
    f = interp1(supp,phi,x,'v5cubic');
    % If the x value was out of the support interp1 returns NaN so we have
    % put zeros everywhere NaN exists.
    f(isnan(f)) = 0;
else % must be mother
    % Interpolate and get function values for input x values.
    f = interp1(supp,psi,x,'v5cubic');
    % If the x value was out of the support interp1 returns NaN so we have
    % put zeros everywhere NaN exists.
    f(isnan(f)) = 0;
end
