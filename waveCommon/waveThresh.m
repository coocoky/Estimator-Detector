%--------------------------------------------------------------------------
% Function:    waveThresh
% Description: Thresholds wavelet coefficients.
%
% Inputs:
%   coeffs            - Nx1 vector of coefficients for the basis functions.
%                       N depends on the number of levels and translations.
%   coeffsIdx         - Lx2 matrix containing the start and stop index
%                       locations of the coeffients for each level in the
%                       coefficient vector.  L is the number of levels.
%                       For example, the set of coefficients for the
%                       starting level can be obtained from the
%                       coefficients vector as:
%                       coeffs(coeffsIdx(1,1):coeffsIdx(1,2),1)
%   startLevel        - starting level for the the father wavelet
%                       (i.e. scaling function).  
%   stopLevel         - last level for mother wavelet scaling.  The start
%                       level is same as the father wavelet's.
%   numSamps          - number of samples used in the density estimation.
% Outputs:
%   coeffs            - Nx1 vector of thresholded coefficients for the 
%                       basis functions.
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

function coeffs = waveThresh(coeffs, coeffsIdx, startLevel, ... 
                             stopLevel, numSamps)
numLevels = length(startLevel:stopLevel);

% Cycle through the levels and threshold.
% Skip the first level b/c we are not thresholding the scaling functions.
levelCoeffs = zeros(numLevels-1,1);
cnt = 0;
% for l = 2 : numLevels
%     cnt = cnt + 1;
%     levelCoeffs(:,cnt) = coeffs(coeffsIdx(l,1):coeffsIdx(l,2));
%     % Alternate: sqrt(l-startLevel)*numSamps^(-.25)
%     levelThresh = sqrt(l-startLevel)/sqrt(numSamps);
%     
% end
if(numLevels>1)
    levelThresh = sqrt((startLevel:stopLevel)-startLevel)/sqrt(numSamps);
else
    levelThresh = 1/sqrt(numSamps);
end


% Compute the proportionality constant by find the one that 
% takes the sum squared coefficients value close to one.
possibleCs = -2:.03:2;
numCs = length(possibleCs);
sqSumTrack = [];

for p=1:numCs
    cnt = 0;
    tempCoeffs = coeffs;
    % apply the level wise thresholds to the coefficients.
    for l = 1:numLevels
        cnt = cnt + 1;
        levelCoeffs = coeffs(coeffsIdx(l+1,1):coeffsIdx(l+1,2));
        levelCoeffs(abs(levelCoeffs)<=possibleCs(p)*levelThresh(cnt))=0;
        tempCoeffs(coeffsIdx(l+1,1):coeffsIdx(l+1,2))=levelCoeffs;
    end
    sqSumTrack = [sqSumTrack norm(tempCoeffs)^2];
end
[val,loc]=max(sqSumTrack);
bestC = possibleCs(loc);

% Threshold with best C.
cnt = 0;
for l = 1:numLevels
    cnt = cnt + 1;
    levelCoeffs = coeffs(coeffsIdx(l+1,1):coeffsIdx(l+1,2));
    levelCoeffs(abs(levelCoeffs)<=bestC*levelThresh(cnt))=0;
    coeffs(coeffsIdx(l+1,1):coeffsIdx(l+1,2))=levelCoeffs;
end

% Renormalize before passing coefficients out.
D = norm(coeffs)^2;
coeffs = coeffs/sqrt(D);
