%--------------------------------------------------------------------------
% Function    : computeRegions
% Description : This function converts the array of domain ranges into an
%               array of domain indices.
% 
% Inputs  :
%     - domainRanges   : array of domain ranges.
%     - discretization : discretization value for the domain.
%     - densitySupport : support of the density function.
% 
% Outputs :
%     - regions : array of domain indices.
%          
% Author(s):
%   - Gedeon Nyengele (gnyengele3@mail.gatech.edu)
%
% Date Published: 6/27/2014
%-------------------------------------------------------------------------- 
function regions = computeRegions( domainRanges,...
                                   discretization,...
                                   densitySupport )
                               
regions = round( (domainRanges - min(densitySupport))./discretization + 1 );

end % function regions = computeRegions( ... )