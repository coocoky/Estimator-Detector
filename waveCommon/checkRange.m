%--------------------------------------------------------------------------
% Function:    checkRange
% Description: Checks to see if the given sample is within the specified
%              range.
%
% Inputs:
%   sample            - Vector of numbers; the sample.
%   sampleDomain      - Maximum and minimum range of the samples. Samples
%                       outside this range will be discarded.
% Outputs:
%   inRange           - 1 if sample value is within range, 0 if not in
%                       range.
% Usage: This function is used to verify whether or not a sample falls in
%        the domain of the density function.
%
% Authors(s):
%   Eddy Ihou(ihouk2002@my.fit.edu)
%   Mark Moyou(mmoyou@my.fit.edu)
%--------------------------------------------------------------------------
function inRange = checkRange(sample, sampleDomain)

inRange = (max(sample) <= max(sampleDomain) && min(sample) >= min(sampleDomain));
if(~inRange)
    warning('Point outside domain. Coefficients not updated.');
end

end % checkRange()
