%--------------------------------------------------------------------------
% Function    : changeDetectorInit
% Description : This function initializes the density change detector.
% 
% Inputs  :
%     NONE
% 
% Outputs :
%     NONE
% 
% Usage   : Use this function to initialize the density change detector.
%           Parameters to Initialise:
%           - CD.NumOfPairs  : sets the number of pairs of k-windows to use
%           - CD.WindowSizes : vector of sizes of each pair of k-windows
%           - CD.Thresholds  : detection thresholds used for the window pairs
%           - CD.Regions     : regions A within which to compare densities.
%                              These regions are defined as intervals.
%           - CD.Dist        : distance measure to use for comparing
%                              densities. This can be set to any of the
%                              distances defined in the DISTANCE CODES
%                              section.
%           
%           This function uses the static (global) CD  (Change Dectector)
%          
% Author(s):
%   - Gedeon Nyengele (gnyengele3@mail.gatech.edu)
%
% Date Published: 6/26/2014
%-------------------------------------------------------------------------- 
function changeDetectorInit()

% Initialize the static (global) Change Detector variable CD.
global CD;
CD = struct();

% CHOOSE THE NUMBER OF WINDOW PAIRS
% -------------------------------------------------------------------------
CD.NumOfPairs = 4;

% CHOOSE THE SIZES OF THE PAIRS OF WINDOWS
% -------------------------------------------------------------------------
CD.WindowSizes =  [200, 400, 600, 800];

% CHOOSE THRESHOLD VALUES FOR EACH WINDOW CHECKING
% -------------------------------------------------------------------------
CD.Thresholds = [0.45, 0.50, 0.40, 0.40];

% CHOOSE REGIONS A OF INTEREST
% -------------------------------------------------------------------------
CD.Regions = [ -2.5 -1.8;...
               -1.0, 1.0 ];
           
% DISTANCE CODES
% -------------------------------------------------------------------------
CD.KSI = 'ksi';
CD.PHI = 'phi';
CD.KS = 'KS';

% CHOOSE A DISTANCE MEASURE
% -------------------------------------------------------------------------
CD.Dist = CD.KSI;

% CONTAINER OF THE TWO DENSITIES TO COMPARE
% -------------------------------------------------------------------------
CD.WindowDensities = cell(CD.NumOfPairs, 2);

% STORES THE SAMPLE COUNTING OF EACH WINDOW PAIR.
% -------------------------------------------------------------------------
CD.WindowContents = ones( CD.NumOfPairs, 2);





end