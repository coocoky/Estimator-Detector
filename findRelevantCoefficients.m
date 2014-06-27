%--------------------------------------------------------------------------
% Function:    findRelevantCoefficients
% Description: Finds the relevant coefficients corresponding to a sample.
%
% Inputs:
%   sample            - Single sample value. 
%   lowerSupp         - Lower value of wavelet support.
%   upperSupp         - Upper value of wavelet support.
%   resolutLevel      - Current resolution level.
%   coefficients      - Vector of coefficients.
%   transVec          - Vector of translate values. 
% 
% Outputs:
%   relevantCoeff     - Relevant coefficients corresponding to sample.
%   relevantKs        - Relevant translates.
%   lowerKIndex       - Lower K index in the translates vector.
%   upperKIndex       - Upper K index in the translates vector.

% Usage: This function is used to find the relevant coefficients in the
%        estimate coefficients function.
%
% Authors(s):
%   Eddy Ihou(ihouk2002@my.fit.edu)
%   Mark Moyou(mmoyou@my.fit.edu)
%--------------------------------------------------------------------------
function [correspCoeff,...
          relevantKs,...
          lowerKIndex,...
          upperKIndex] = findRelevantCoefficients(sample,...
                                                  lowerSupp,...
                                                  upperSupp,...
                                                  resolutLevel,...
                                                  coefficients,...
                                                  transVec)

    % Find the translate bounds that correspond to the current sample.
    kLowerBound = ceil(2^resolutLevel*sample - upperSupp);
    kUpperBound = floor(2^resolutLevel*sample - lowerSupp);
    
    % Find the translates in the translate vector that correspond to the    
    % sample.
    lowerKIndex = find(transVec == kLowerBound);
    upperKIndex = find(transVec == kUpperBound);
    relevantKs = (transVec(lowerKIndex) : transVec(upperKIndex));
    
    % Find the corresponding coefficients. 
    correspCoeff = coefficients(lowerKIndex : upperKIndex);
    