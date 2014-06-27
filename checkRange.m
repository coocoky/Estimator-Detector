
function inRange = checkRange(sample, domain)

inRange = (max(sample) <= max(domain) && min(sample) >= min(domain));


end % checkRange()
