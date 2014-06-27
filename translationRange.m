
function transRange = translationRange(level)

global DensityVars;

wSupport      = DensityVars.WaveSupport;
transRange    = zeros(1,2);
transRange(1) = floor((2^level)*DensityVars.DensityDomain(1)-wSupport(2));
transRange(2) = ceil((2^level)*DensityVars.DensityDomain(2)-wSupport(1));

