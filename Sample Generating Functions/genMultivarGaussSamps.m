%--------------------------------------------------------------------------
% Function:    genMultivarGaussSamps
% Description: Generates multivariate Gaussian samples for a given mean
%              and covariance matrix.
%
% Inputs:
%   numSamps   - number of samples you want to draw from the Gaussian, 
%                i.e. x will be elements in R^n.
%   mu         - nX1 vector of means where n is the dimension of the
%                space, i.e. mu element in R^n.
%   covariance - nXn covariance matrix for the mulitvariate Gaussian, 
%                where n is the dimension of the space.
%
% Outputs:
%   samps      - numSampsXn matrix of samples
%
% Usage:
%   Create 10,000 samples with mean centered at (0,1) and with 
%   x-axis variance equal 1 and y-axis variance equal 10.
%   >> covariance = [1 0;0 10];
%   >> mu = [0 1];
%   >> numSamps = 10e3;
%   >> samps = genMultivarGaussSamps(numSamps, mu, covariance);
%
% Authors(s):
%   Adrian M. Peter
%--------------------------------------------------------------------------
function samps = genMultivarGaussSamps(numSamps, mu, covariance,varargin)

mu = mu(:)';

n = size(mu,2);

% Create n dimensional, mean zero samples.
samps = randn(numSamps,n);
if(numSamps>1) %Doens't always come out w/ zero mean so will remove.
    tempMu = mean(samps);
    samps = samps - repmat(tempMu,numSamps,1);
end
% Perform Cholesky decomp of the covariance matrix.
% Any positive semidefinite matrix A can be factored in the form A=LL' 
% for some real square matrix L, which we may think of as a matrix 
% square root of A.  One factorization to perform this known as the 
% Cholesky factorization. Any positive semidefinite matrix A has a 
% factorization of the form A=LL' where L is a lower triangular matrix.
% This is exactly what we want in order for our samples to have the right
% covariance.  If the A is positive definite, then the call to
% chol returns p=0.  But if A is not positive definite, then p is a 
% positive integer, in which case we can use sqrtm.
[R,p]=chol(covariance);
if p>0
    samps=samps*sqrtm(covariance);
else
    samps=samps*R;
end

% Add the mean
samps = samps + repmat(mu,numSamps,1);

%if(isempty(varargin(1)))
    plotFlag = 1;
%else
%    plotFlag = cell2mat(varargin(1));
%end
if(plotFlag)
    if(size(samps,2)==2)
        plot(samps(:,1),samps(:,2),'ro');
    else
        plot(samps,'ro');
    end
end
