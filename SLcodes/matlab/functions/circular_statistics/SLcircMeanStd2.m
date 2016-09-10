
% SLcircMeanStd2.m
%
%  author: Steeve Laquitaine
% purpose: Calculate the circular statistics of the data
%
%   usage:  
%
%        data = SLcircMeanStd2([0,1; 1,0])
%        data = SLcircMeanStd2([90;315],'polar')
%        data = SLcircMeanStd2([45;45],'polar')
%important: when y are polars, y must be a column vector


function data = SLcircMeanStd2(y,varargin)

%case y is polar (deg)
%convert to cartesian
if size(y,2)==1 || sum(strcmp(varargin,'polar'))==1
    radius = ones(size(y,2));
    coord = polar2cartesian(y,radius);
elseif size(y,2)~=1 || sum(strcmp(varargin,'cartesian'))==1
    %cartesian
    coord=y;
end

%register the coordinates of the input directions
data.coord.all = coord;

%convert from cartesian coordinates to angles (in degree)
data.deg.all = getangle(coord(:,1),coord(:,2));

%cartesian mean and mean (deg)
data.coord.mean = nanmean(coord,1);
data.deg.mean = getangle(data.coord.mean(:,1),data.coord.mean(:,2));

%calculate the std to the mean direction est (in degree);!!! could be a
%subfunction itself.
%Apply the rule of thumb that follows. It seems to work fine intuitively. 
%It would be nice to fine a cleaner way to calculate the std.
%initialize the 'sample' and 'mean' variables used to calculate the std
data.num =  numel(data.deg.all);
data.deg.allforstd = data.deg.all;
data.deg.meanforstd = repmat(data.deg.mean,data.num,1);

%if the resulting mean direction is between 0 and 180.
if data.deg.mean+180 <= 360
    
    %if estimation is>=mean direction+180
    for i = 1 : data.num
        if data.deg.all(i) >= data.deg.mean+180
            data.deg.allforstd(i) = data.deg.all(i) - 360;
        end
    end
    
else
    
    %if the resulting mean direction is between 180 and 360.
    %if the estimated direction sampled is <=the mean direction - 180
    %This operation permits to linearized distance error for var,std and
    %sem
    %360 is substracted from angle > 180
    for i = 1 : data.num
        if data.deg.all(i) <= data.deg.mean - 180
            data.deg.meanforstd(i) = data.deg.mean - 360;
        end
    end
    
end

%variance, std and sem in term of distance to the mean
%Distance error is linearized
%360 is substracted from angle > 180
%0 is substracted from angle < 180
distanceError = (data.deg.allforstd - data.deg.meanforstd);
data.deg.var = nanmean(distanceError.^2,1);
data.deg.std = sqrt(data.deg.var);
data.deg.sem = data.deg.std/sqrt(data.num);


