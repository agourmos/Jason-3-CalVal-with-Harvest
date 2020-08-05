function [distances] = minDist(lat,lon)
% This function finds the distance between 2 points on a map given a
% longitude and latitude
    LON = -120.671667*(pi/180);
    LAT = 34.468333*(pi/180);
    R = 6371; %Km
    distances = zeros(1,length(lat));
    for i = 1:length(lat)
        conValLat = lat(i)*(pi/180);
        conValLon = lon(i)*(pi/180);
        a = (sin((conValLat-LAT)/2)^2)+(cos(LAT)*cos(conValLat)...
            *(sin(((conValLon-LON)/2)^2)));
        c = 2 * atan2(sqrt(a), sqrt(1-a));
        distances(i) = R*c;
    end
end

