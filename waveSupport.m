%--------------------------------------------------------------------------
% Function:    waveSupport
% Description: Returns the support the given wavelet.
%
% Inputs:
%   wName             - 3 to 4 charater code name of wavelet for density 
%                       approximation.  Currently supports:
%                       'db1', 'db2', 'db3', 'db4', 'db5'
%                       'db6', 'db7', 'db7', 'db8', 'db9','db10'
%                       'meyr'
%                       Note: db1 is Haar.
% Outputs:
%   support           - 1x2 vector of the wavlet support.
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

function support = waveSupport(wName)
switch(lower(wName))
    case {'db1';'db2';'db3';'db4';'db5';'db6';'db7';'db8';'db9';'db10'}
        bLoc = strfind(wName,'b');
        N = str2num(wName(bLoc+1:end));
        support = [0 2*N-1];
    case 'dmey'
        support = [0 101];
    case {'sym1';'sym2';'sym3';'sym4';'sym5';'sym6';'sym7';'sym8';'sym9';'sym10'}
        bLoc = strfind(wName,'m');
        N = str2num(wName(bLoc+1:end));
        support = [0 2*N-1];
    case {'coif1';'coif2';'coif3';'coif4';'coif5';'coif6';'coif7';'coif8';'coif9';'coif10'}
        bLoc = strfind(wName,'f');
        N = str2num(wName(bLoc+1:end));
        support = [0 6*N-1];
    otherwise
        error('Invalid wavelet type! Please do help on ''wavesupport'' for valid list!');
end
