
% SLcircMeanStd.m
%
%  author: Steeve Laquitaine
% purpose: Calculate the circular statistics of the data
%
%   usage:  
%
%        data = SLcircMeanStd([0,1; 1,0])
%        data = SLcircMeanStd([0; 90])
%        data = SLcircMeanStd([0; 90],'polar')
%important: when y are polars, y must be a column vector


function data = SLcircMeanStd(y,type)

%check type
if nargin<2
    
    fprintf('\n ------------(Missing Argument)------------ \n')
    fprintf('(SLcircMeanStd) Please indicate data type  \n')
    fprintf('(SLcircMeanStd) You choices are :          \n')
    fprintf('(SLcircMeanStd) - "polar"                  \n')
    fprintf('(SLcircMeanStd) - "cartesian"              \n')
    dbstack
    keyboard
    
end

%case y is polar (deg)
%convert or get cartesian data
if size(y,2)==1 || sum(strcmp(type,'polar'))==1    
    radius = ones(size(y,2));
    coord = SLpolar2cartesian(y,radius,'polar');    
elseif size(y,2)~=1 || sum(strcmp(type,'cartesian'))==1
    coord = y;
end

%directions cartesians
data.coord.all = coord;

%convert cartesian to angles (in degree)
data.deg.all = SLgetangle(coord(:,1),coord(:,2));

%cartesian mean and mean (deg)
%the mean will be 0,0 if the two angles are opposite
%and the mean angle is undefined : nan
data.coord.mean = nanmean(coord,1);
data.deg.mean = SLgetangle(data.coord.mean(:,1),data.coord.mean(:,2));

%tag special cases when means are undefined
%because angles cancel each other : 
%mean coordinates are 0,0 and
%mean angle is null            
if data.coord.mean(1)==0 && data.coord.mean(2)==0
    data.deg.isMeanNull = 1;
else
    data.deg.isMeanNull = 0;
end
    
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
    for i = 1 : data.num
        if data.deg.all(i) <= data.deg.mean - 180
            data.deg.meanforstd(i) = data.deg.mean - 360;
        end
    end
    
end

%variance, std and sem
%variance, std and sem in term of distance to the mean
%Distance error is linearized
data.deg.var = nanmean((data.deg.allforstd - data.deg.meanforstd).^2,1);
data.deg.std = sqrt(data.deg.var);
data.deg.sem = data.deg.std/sqrt(data.num);


