
% vectorStat.m
% 
% author: steeve laquitaine
%    date: 131128
% 
%   usage: data=vectorStat([0 1; 0.5 0.2; 0.1 1])
%           "coord" are cartesian coordinates of the vector;
%           origin of the vector at 0 and end point at "coord".
% 
% purpose: calculate statistics of directional data (vectors)


function data=vectorStat(coord)

%convert from cartesian to angles (in degree)
data.coord.all=coord;
if isempty(coord)
    coord=[nan nan];
end
data.deg.all=SLgetangle(coord(:,1),coord(:,2));

%mean vector's coordinates
data.coord.mean=nanmean(coord,1);

%mean angle in degree
data.deg.mean=SLgetangle(data.coord.mean(:,1),data.coord.mean(:,2));

%mean angle's variance, std & sem
data.num=numel(data.deg.all);
data.deg.allforstd=data.deg.all;
data.deg.meanforstd=repmat(data.deg.mean,data.num,1);

%if the resulting mean direction is between 0 and 180.
if data.deg.mean+180<=360
    %if estimation is>=mean direction+180
    for i=1:data.num
        if data.deg.all(i)>=data.deg.mean+180
            data.deg.allforstd(i)=data.deg.all(i)-360;
        end
    end
    %if the resulting mean direction is between 180 and 360.
else
    %if the estimated direction sampled is <=the mean direction - 180
    for i=1:data.num
        if data.deg.all(i)<=data.deg.mean-180
            data.deg.meanforstd(i)=data.deg.mean-360;
        end
    end
end
data.deg.var=nanmean((data.deg.allforstd-data.deg.meanforstd).^2,1);
data.deg.std=sqrt(data.deg.var);
data.deg.sem=data.deg.std/sqrt(data.num);