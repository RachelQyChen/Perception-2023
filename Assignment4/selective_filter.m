%Side-selective filter
function [evenLeft, oddLeft, evenRight, oddRight] = selective_filter(oddFast, oddSlow, evenFast, evenSlow)  
    evenLeft = oddFast + evenSlow;
    oddLeft = -oddSlow + evenFast;
    evenRight = -oddFast + evenSlow;
    oddRight = oddSlow + evenFast; 
end