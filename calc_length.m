function y = calc_length(tspan,y,maxy)
%this function checks if the bee reached maximum temperature before the end
%of tspan 

if length(y)<length(tspan)
    y((length(y)+1):length(tspan))=maxy;   %if it got cut off with T_th>maxy, fill in the rest of y
end

