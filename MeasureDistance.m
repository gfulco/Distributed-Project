function [ measure ] = MeasureDistance(self,other, maxrange )
% This function checks if there is an agent near the another
% If the other agent is in the range of sensors then the relative position
% is calculated and the parameters for the filter are returned

relx = (other(1)-self(1));
rely = (other(2)-self(2));
measure = (sqrt(relx^2+rely^2));

if isnan(maxrange)
    measure = measure + (2*rand()-1)*0.2;
elseif(measure <= maxrange)
    measure = measure + (2*rand()-1)*0.02;
else
    measure = NaN;
end
end

