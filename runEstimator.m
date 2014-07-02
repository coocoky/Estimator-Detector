%--------------------------------------------------------------------------
% Script      : runEstimator
% Description : This script runs the density estimator. It also runs the
%               density change detector if this functionality is enabled.
%          
% Author(s):
%   - Gedeon Nyengele (gnyengele3@mail.gatech.edu)
%
% Date Published: 6/26/2014
%-------------------------------------------------------------------------- 
clear all; clc; close all

addpath('Generated Samples');
addpath('Sample Generating Functions');
addpath('waveCommon');
addpath('TestDataSets');

global DensityVars CD;

% Initialize density estimator variables.
estimatorInit();

% Initialize translates.
getTranslates();

% Initialize coefficients.
initializeCoefficients();

% Run online density estimator.
densityDomainPlot = (DensityVars.DensityDomain(1) : DensityVars.Discretization : DensityVars.DensityDomain(2));
sampleLength = length(DensityVars.Sample);
sampBforWindow = DensityVars.DensityDomain(1) - 1;

% Initialize the density change detector if enabled.
if( DensityVars.ChangeDetectorFlag == 1)
    changeDetectorInit();
end % if( DensityVars.ChangeDetectorFlag == 1)


for i = 1 : sampleLength
    % ---------------------------------------------------------------------
    %                  GARCIA-TREVINO WINDOW-METHOD
    % ---------------------------------------------------------------------
    
    % Update sample leaving window.
    if (i > DensityVars.WindowSize)
        sampBforWindow = DensityVars.Sample((i - DensityVars.WindowSize) + 1);
    end % if (i > DensityVars.WindowSize)
    
    % Estimating DensityVars.Coefficients.
    estimateCoefficients(DensityVars.Sample(i),i, sampBforWindow); 
    
    % Showing plots of density. showPlots = 1 --- on.
    if (DensityVars.ShowPlot == 1 && mod(i, DensityVars.PlotUpdateFreq)== 0)
        
            % Obtaining the updated density for each sample.
            density = calculateDensity();

            % Clear figure.
            clf;

            plot(densityDomainPlot, density, 'r', 'LineWidth', 1.5);
            title(['Samples: ', num2str(i)]);
            xlabel('X');
            ylabel('Density');
            
            % Setting y limit to change for demoType ==
            % 'strongSkewUniToClaw'.
            if (strcmp(DensityVars.DemoType,'strongSkewUniToClaw'))
                if (i < DensityVars.NumSampsDistr1)
                    ylim([0 1.5]);
                end
            else
                ylim([0 0.75]);
            end
            box on;
            drawnow();
             
    end %   if (DensityVars.ShowPlots == 1 && mod(i, DensityVars.PlotUpdateFreq)== 0)
    
    
    
    % ---------------------------------------------------------------------
    %               KIFER WINDOWS FOR CHANGE DETECTION
    % ---------------------------------------------------------------------
    
    if( DensityVars.ChangeDetectorFlag == 1)
        
        for h = 1:CD.NumOfPairs
            
            % Fill the static window of the pair if it's not full.
            if( CD.WindowContents(h, 1) <= CD.WindowSizes(h) )
                CD.WindowContents(h, 1) = CD.WindowContents(h, 1) + 1;
                
                % If the static window is full, record its density.
                if( CD.WindowContents(h, 1) == CD.WindowSizes(h) )
                    CD.WindowDensities{h, 1} = calculateDensity();
                end % if( CD.WindowContents(h, 1) == CD.WindowSizes(h) )
                
            
            % Fill, the moving window if not full.
            elseif( CD.WindowContents(h, 2) < CD.WindowSizes(h) )
               
                CD.WindowContents(h, 2) = CD.WindowContents(h, 2) + 1; 
            
            % If the two windows are full, 
            else
                % Record the density of the sliding window.
                CD.WindowDensities{h, 2} = calculateDensity();
                
                % Try to detect a change.
                regions = computeRegions( CD.Regions, DensityVars.Discretization, DensityVars.DensityDomain );
                [changed, dist] = detectChange(CD.WindowDensities{h,1}, CD.WindowDensities{h,2}, CD.Thresholds(h), CD.Dist, regions);
                
                % If density changed
                if( changed == 1 )
                    
                    % Report the change.
                    reportChange(i, h, dist);
                    
                    % Clear windows that detected the change.
                    % CD.WindowContents(h,:) = [1, 1];
                    CD.WindowContents = ones(CD.NumOfPairs, 2);
                end % if( changed == 1 )
            end % if( CD.WindowContents(h, 1) < CD.WindowSizes )
            
            
        end % for h = 1:CD.NumOfPairs
        
    end % if( DensityVars.ChangeDetectorFlag == 1)
    
end % i = 1 : sampleLength





