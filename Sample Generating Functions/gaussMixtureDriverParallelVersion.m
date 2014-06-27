%--------------------------------------------------------------------------
% Script:    gaussMixtureDriver
% Description: Creates different types of mixture models.
%              Models come from:
%               @article{Marron92,	
%                   title= "Exact mean integrated squared error",
%                   author= "S. J. Marron and M. P. Wand",
%                   journal= "The Annals of Statistics",
%                   volume="20",
%                   number="2",
%                   pages="712--736"
%                   year="1992",
%               }
%
% Authors(s):
%   Adrian M. Peter
%--------------------------------------------------------------------------
close all; clc;
clear all;

% GMM parameters.
% Supported names:
% gaussian, skewedUni, strongSkewedUni, kurtoticUni, outlier, bi,
% separatedBi, skewedBi, tri, claw, dblClaw, asymmClaw, asymmDblClaw
gmmName    = 'claw';
numSamps   = 10e6;
genSamples = 1;
saveFileName=[gmmName num2str(numSamps)];

% Opening matlabpool
runParallel = 1;
% Set number of workers
noWorkers = 8;
isOpen = matlabpool('size') > 0;
isOpenCorr = matlabpool('size') == noWorkers;
if runParallel && ~isOpenCorr,
    matlabpool close force
    matlabpool(noWorkers)
end
if ~runParallel && isOpen,
    matlabpool close force
end
tic;
% Based on name setup the GMM priors, means, variances.
switch(gmmName)
    case 'gaussian'
        n = 1;
        m = @(x) 0;
        v = @(x) 1;
        p = @(x) 1;
        gmmFullName = 'Gaussian';
    case 'skewedUni'
        n = 3;
        m = @(x) (x==1)*0 + (x==2)*.5 + (x==3)*13/12;
        v = @(x) (x==1)*1 + (x==2)*(2/3)^2 + (x==3)*(5/9)^2;
        p = @(x) [1/5 1/5 3/5];
        gmmFullName = 'Skewed Unimodal';
    case 'strongSkewedUni'
        n = 8;
        m = @(x) 3*((2/3)^(x-1)-1);
        v = @(x) (2/3)^(2*(x-1));
        p = @(x) [1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8];
        gmmFullName = 'Strongly Skewed Unimodal';
    case 'kurtoticUni'
        n = 2;
        m = @(x) 0;
        v = @(x) (x==1)*1 + (x==2)*(1/10)^2;
        p = @(x) [2/3 1/3];
        gmmFullName = 'Kurtotic Unimodal';
    case 'outlier'
        n = 2;
        m = @(x) 0;
        v = @(x) (x==1)*1 + (x==2)*(1/10)^2;
        p = @(x) [1/10 9/10];
        gmmFullName = 'Outlier';
    case 'bi'
        n = 2;
        m = @(x) (x==1)*-1 + (x==2)*1;
        v = @(x) (x==1)*(2/3)^2 + (x==2)*(2/3)^2;
        p = @(x) [.5 .5];
        gmmFullName = 'Bimodal';
    case 'separatedBi'
        n = 2;
        m = @(x) (x==1)*-3/2 + (x==2)*3/2;
        v = @(x) (x==1)*(1/2)^2 + (x==2)*(1/2)^2;
        p = @(x) [.5 .5];
        gmmFullName = 'Separated Bimodal';
    case 'skewedBi'
        n = 2;
        m = @(x) (x==1)*0 + (x==2)*3/2;
        v = @(x) (x==1)*1 + (x==2)*(1/3)^2;
        p = @(x) [.75 .25];
        gmmFullName = 'Skewed Bimodal';
    case 'tri'
        n = 3;
        m = @(x) (x==1)*-(6/5) + (x==2)*(6/5) + (x==3)*0;
        v = @(x) (x==1)*(3/5)^2 + (x==2)*(3/5)^2 + (x==3)*(1/4)^2;
        p = @(x) [9/20 9/20 1/10];
        gmmFullName = 'Trimodal';
    case 'claw'
        n = 6;
        m = @(x) ((x-1)==0)*0 + ((x-1)>0)*(((x-2)/2)-1);
        v = @(x) ((x-1)==0)*1 + ((x-1)>0)*(1/10)^2;
        p = @(x) [.5 1/10 1/10 1/10 1/10 1/10];
        gmmFullName = 'Claw';
    case 'dblClaw'
        n = 9;
        m = @(x) (x==1)*-1 + (x==2)*1 + (x==3)*-3/2 + (x==4)*-1 ...
                +(x==5)*-.5 + (x==6)*0 + (x==7)*.5 + (x==8)*1 + (x==9)*1.5;
        v = @(x) (x==1)*(2/3)^2 + (x==2)*(2/3)^2 + (x>2)*(1/100)^2
        p = @(x) [49/100 49/100 1/350 1/350 1/350 1/350 1/350 1/350 1/350];
        gmmFullName = 'Double Claw';
    case 'asymmClaw'
        n = 6;
        m = @(x) (x==1)*0 + (x==2)*-1.5 + (x==3)*-.5 ...
                +(x==4)*.5 + (x==5)*1.5 + (x==6)*2.5;
        v = @(x) (x==1)*1 + (x==2)*(2^2/10)^2 + (x==3)*(2/10)^2 ...
                +(x==4)*(1/10)^2 + (x==5)*(2^-1/10)^2 + (x==6)*(2^-2/10)^2;
        p = @(x) [.5 (2^3)/31 (2^2)/31 (2/31) (1/31) (2^-1)/31];
        gmmFullName = 'Asymmetric Claw';
    case 'asymmDblClaw'
        n = 8;
        m = @(x) (x==1)*-1 + (x==2)*1 + (x==3)*-.5 + (x==4)*-1 ...
                +(x==5)*-1.5 + (x==6)*.5 + (x==7)*1 + (x==8)*1.5;
        v = @(x) (x==1)*(2/3)^2 + (x==2)*(2/3)^2 + ((x>2)&(x<=5))*(1/100)^2 ...
                +((x>5)&(x<=8))*(7/100)^2;
        p = @(x) [46/100 46/100 1/300 1/300 1/300 7/300 7/300 7/300];
        gmmFullName = 'Asymmetric Double Claw';
    case 'smoothComb'
    case 'discreteComb'
    otherwise
        error('Invalid GMM Type!');
end

if(genSamples)
    % Generate the 
    % Set the state to make results reproducible
    randn('state',316);
    rand('state', 316);
    priorCDF = [0 cumsum(p(n))];
    samps = zeros(numSamps,1);
    parfor s = 1 : numSamps
        unifRv = rand;
        for comp = 2 : n + 1
            if ((unifRv >= priorCDF(comp - 1)) & (unifRv < priorCDF(comp)))
                samps(s) = genMultivarGaussSamps(1, m(comp-1), v(comp-1),0);
            end
        end
    end
end
% Create analytic plot of GMM.
syms x;
gmm = 0;
gmmTruth = 0;
xx = [-3.5:.01:3.5];
for comp = 1 : n
    pp = p(n);
    gmm = gmm + pp(comp)*(1/sqrt(2*pi*v(comp)))*exp(-.5*((x-m(comp))^2)/v(comp)); 
    gmmTruth = gmmTruth + pp(comp)*(1/sqrt(2*pi*v(comp)))*exp(-.5*((xx-m(comp)).^2)/v(comp)); 
end
%figure;ezplot(gmm,[-3.5 3.5]);
figure;plot(xx,gmmTruth);

title(gmmFullName,'FontSize',14);
if(genSamples)
    figure;plot(samps,'ro');
    title('Samples','FontSize',14);
    domain = xx;
    save(['Z:\Users\mmoyou\BensWaveletDensityEstimators\Density Estimator Ammenities\Generated Samples\' saveFileName],'samps','gmmTruth','domain');
    %save(['../../MyTutorials/WaveletDensityEstimation/dynamical/1D/' gmmName],'samps','gmmTruth');
end
toc;