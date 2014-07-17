Estimator-Detector
==================

Density estimator with density-change detection using wavelets.

##How to use the density change detector
This 1D version of the density change detector comes bundled with the Matlab 1D density estimator.
The changed detector is enabled by default. 

To disable the change detector:
- open the changeDetectorInit.m file
- Inside the initialization function, set CD.ENABLE to 0.
The change detector is disabled when the value of CD.ENABLE is set to 0 and it is enabled when the value of CD.ENABLE is set to 1.


The current implementation of the change detector uses pairs of sliding windows and continuously compares the enclosed approximated density function by the two windows in the pair.
If the two density functions differ significantly, the detector flags a change.
It is important to note that the measure of similarity or dissimilarity between the density functions approximated by the two given windows in the pair depends on the distance measure used.
With the current version of the change detector, we provide 3 distance measures: ksi, phi, and KS distances.

The number of pairs of windows to used can be set by changing the CD.NumOfPairs field within the change detector's initialization function.
After defining the number of pairs of windows to use, it is essential to also set the sizes for the pairs of windows. The windows in the same pair get the same size.
So the CD.WindowSizes field defines the window sizes for each pair of windows.
Equally important is the need to set the threshold distances for the change detector's flag.
Please set a threshold distance per pair of windows by changing the CD.Thresholds field.

If one is interested in detecting changes in particular regions of the domain of the density function, these regions can be defines in the CD.Regions field.
for each region, please specify the two endpoints for the region as a row entry in the CD.Regions field.
