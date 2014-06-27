%--------------------------------------------------------------------------
% Function:    translationRange
% Description: Calculates the span of translations for the current level.
%              The range of translations is calculated such that one full
%              wavelet support starts before the minimum sample support and
%              one full wavelet support comes after the maximum sample
%              support value.  
%
% Inputs:
%   sampleSupport     - 2x1 vector containing the sample support.
%   wName             - 3 to 4 charater code name of wavelet for density 
%                       approximation.
%   level             - integer value for the wavelet scale level.
%                       Scale is calculated as 2^level.
%   
% Outputs:
%   transRange        - 1x2 vector containing the start and stop
%                       translation values.
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
function transRange = translationRange(sampleSupport, wName, level)
wSupport      = waveSupport(wName);
transRange    = zeros(1,2);
transRange(1) = floor((2^level)*sampleSupport(1)-wSupport(2));
transRange(2) = ceil((2^level)*sampleSupport(2)-wSupport(1));

