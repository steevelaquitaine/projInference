% Calculate an vector average
function [data] = vectorAverage(coord)
% Calculate the cartesian coordinates of the mean estimated direction.
data.coord.mean = nanmean(coord,1); % cartesian coordinate of the mean direction est

% Calculate the mean direction estimate (in degree)
x = data.coord.mean(:,1);
y = data.coord.mean(:,2);

% radian
data.radian.mean = atan(y./x); %theta = arctan(opposite/adjacent);
% degrees
data.degree.mean = data.radian.mean.* 180/pi;
