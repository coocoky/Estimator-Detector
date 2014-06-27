%--------------------------------------------------------------------------
% Function    : estimatorInit
% Description : This function initializes the density estimator variables.
% 
% Inputs  :
%     NONE
% 
% Outputs :
%     NONE
% 
% Usage   : The function builds up the global (static) variable 'DensityVars'
%           with all the different fields with each field being a variable of
%           of the density estimator.
%           Those fields are: 
%           - Samples             : contains the sample dataset.
%           - DensityDomain       : domain of the density function
%           - Wavelet             : 3-4 character codename for the type of
%                                   wavelet used.
%           - StartLevel          : initial resolution of the multiresolution
%                                   analysis.
%           - StopLevel           : final resolution of the multiresolution
%                                   analysis.
%           - WaveletFlag         : flag to enable/disablle mother wavelet.
%                                   0 -> mother wavelet disabled.
%                                   1 -> mother wavelet enabled.
%           - WindowSize          : size (length) of the window used for the
%                                   Garcia-Trevino's window method
%           - Discretization      : discretization value for the density do-
%                                   main.
%           - AgingMethod         : method used for the aging of old samples.
%                                   0 -> No Aging.
%                                   1 -> Window Method.
%                                   2 -> Caudle Aging.
%           - PlotUpdateFreq      : frequency at which to plot the density.
%           - ShowPlot            : flag to enable/disable density plotting.
%                                   0 -> Plotting disabled.
%                                   1 -> Plotting enabled.
%           - PhiInterpolant      : gridded interpolant class used to speed up
%                                   interpolation of scaling(phi) function.
%           - PsiInterpolant      : gridded interpolant class used to speed up
%                                   interpolation of wavelet(psi) function.
%          
% Author(s):
%   - Gedeon Nyengele (gnyengele3@mail.gatech.edu)
%
% Date Published: 6/26/2014
%--------------------------------------------------------------------------
function estimatorInit()

% Initialize the static (global) DensityVars.
global DensityVars;
DensityVars = struct();

DensityVars.AgingFlag          = 2;
DensityVars.WindowSize         = 600;
DensityVars.Wavelet            = 'db4';
DensityVars.StartLevel         = 0; 
DensityVars.StopLevel          = 2; 
DensityVars.WaveletFlag        = 1;
DensityVars.Discretization     = 0.1;
DensityVars.DensityDomain      = [-3.5 3.5];
DensityVars.ShowPlot           = 1; 
DensityVars.PlotUpdateFreq     = 25;
DensityVars.ChangeDetectorFlag = 1;
DensityVars.WaveSupport        = waveSupport( DensityVars.Wavelet );

load(['waveCommon/', DensityVars.Wavelet, 'Tables.mat'], 'supp', 'phi', 'psi');
DensityVars.PhiInterpolant = griddedInterpolant(supp, phi, 'cubic');
DensityVars.PsiInterpolant = griddedInterpolant(supp, psi, 'cubic');

DensityVars.DemoType = 'strongSkewUniToClaw';

switch (DensityVars.DemoType)
    case 'claw1000'    
        load('claw1000.mat', 'samps');
        % load('claw10000000.mat', 'samps');
        DensityVars.Sample = samps;
        DensityVars.NumSampsDistr1 = 0;
    case 'gauss2Bimodal2000'
        % Gaussian to Bimodal.(2000 each)
        load('gauss2kToBimodal2k.mat')
        DensityVars.Sample = samps;
        DensityVars.NumSampsDistr1 = sampsPerDistr(1);
    case 'claw10000'
        load('claw10000.mat', 'samps');
        DensityVars.Sample = samps;
        DensityVars.NumSampsDistr1 = 0;
    case 'gauss2Claw2000'
        % Gaussian to claw (2000 each)
        load('gaussian2kToClaw2k.mat')
        DensityVars.Sample = samps;
        DensityVars.NumSampsDistr1 = sampsPerDistr(1);
    case 'strongSkewUniToClaw'
        % strongSkewUni(3000) to claw(4000)
        load('strongSkewUni3kToClaw4k2.mat')
        DensityVars.Sample = samps;
        DensityVars.NumSampsDistr1 = sampsPerDistr(1);
     case 'gaussianToSeperatedBi'
        % gaussian(2000) to seperatedBi(3s000)
        load('seperatedBi3kToGaussian2k.mat')
        DensityVars.Sample = samps;
        DensityVars.NumSampsDistr1 = sampsPerDistr(1);   
    otherwise
        fileToLoad = horzcat(DensityVars.DemoType, '.mat');
        if( exist(fileToLoad, 'file') )
            load( fileToLoad );
            if( exist('samps', 'var') )
                DensityVars.Sample = samps;
            else
                error('Could not find variable samps in the file.');
            end % if( exist('samps', 'var') )
        else 
            error('Invalid data stream!');
        end % if( exist(fileToLoad, 'file') )
end % switch (DensityVars.DemoType)

end % function estimatorInit()
