%--------------------------------------------------------------------------
% Function:    detectChange
% Description: Determines whether or not a change has occurred between
%               the reference and the current density.
%
% Inputs:
%   refDensity        - Vector, the reference density to be compared to.
%
%   curDensity        - Vector, the current density distribution.
%
%   threshold         - Double, the level of the statistic above which
%                       a change should be reported.
%   distanceMeasure   - String, the type of distance measure to use between
%                       the two densities.
%   AIntervals        - N x 2 matrix, the intervals over which to calculate
%                       the distance between densities. Given as indices
%                       into the density vectors.
%   
%
% Outputs:
%   changed            - Boolean, whether or not a change was detected.
% Authors(s):
%   Dan Weinand(daniel.weinand@pomona.edu)
%   Gedeon Nyengele
%
% Date: June 26th, 2014
%--------------------------------------------------------------------------
function [ changed, maxDistance ] = detectChange(refDensity,...
                                                  curDensity,...
                                                  threshold,...
                                                  distanceMeasure,...
                                                  AIntervals)
                                
    
    % Determine the distance
    maxDistance = 0;
    
    for A = AIntervals'
        % Reference and current densities under A
        refADensity = sum(refDensity(A(1):A(2)));
        curADensity = sum(curDensity(A(1):A(2)));
        
        switch distanceMeasure
            case 'KS'
                distance = 2*abs(refADensity - curADensity);
            case 'ksi'
                top = abs(refADensity - curADensity);
                bot = sqrt(min((refADensity + curADensity)/2, 1 - (refADensity + curADensity)/2));
                distance = top/bot;
            case 'phi'
                top = abs(refADensity - curADensity);
                bot = sqrt((refADensity + curADensity)/2 * (1 - (refADensity + curADensity)/2));
                distance = top/bot;
            otherwise
                error('Invalid distance measure');
        end
        
        maxDistance = max(maxDistance, distance);
    end
    
    % The density has changed if the distance is greater than the threshold
    changed = maxDistance > threshold;
end