%SLstatcircular.m
%
%   author: steeve laquitaine
%     date: 140107
%  purpose: calculate circular data mean, variance and std
%
%    usage: 
%
%       data = SLstatcircular([0 1;1 0])
 
function data = SLstatcircular(coord)


%register the coordinates of the input dirs
data.coord.all = coord;

%convert from cartesian coordinates to angles (in degree)
data.deg.all = SLgetangle(coord(:,1), coord(:,2));

%calculate the cartesian coordinates of the mean estimated dir.
data.coord.mean(:,1) = nanmean(coord(:,1)); %cartesian coordinate of the mean dir est
data.coord.mean(:,2) = nanmean(coord(:,2));

%calculate the mean dir estimate (in degree)
data.deg.mean = SLgetangle(data.coord.mean(:,1),data.coord.mean(:,2)); %(in degree)

%calculate the std to the mean dir est (in degree)
%Apply the rule of thumb that follows. It seems to work fine intuitively.
%It would be nice to fine a cleaner way to calculate the std.
%initialize the 'sample' and 'mean' variables used to calculate the std
%count the number of data
data.num=numel(data.deg.all); %sample size
%collect data for the calculation of the std
data.deg.allforstd=data.deg.all;
%collect the mean of the data for the calculation of the std
data.deg.meanforstd=repmat(data.deg.mean,data.num,1);

%if the resulting mean dir is between 0 and 180.
%if data.deg.mean + 180 <= 360
if data.deg.mean <= 180
    %sample each est dir
    for i=1:data.num
        %if the est dir sampled is >=mean dir + 180
        if data.deg.all(i) >=data.deg.mean + 180
            data.deg.allforstd(i)=data.deg.all(i) - 360;
        end
    end
    %if the resulting mean dir is between 180 and 360.
else
    %sample each est dir
    for i=1:data.num
        %if the est dir sampled is <=the mean dir - 180
        if data.deg.all(i) <=data.deg.mean-180
            data.deg.meanforstd(i) = data.deg.mean - 360;
        end
    end
end

%now calculate the variance of the estimated dir.
data.deg.var = nanmean((data.deg.allforstd - data.deg.meanforstd).^2,1);

%and now calculate the std
data.deg.std = sqrt(data.deg.var);

%and now calculate the sem
data.deg.sem = data.deg.std/sqrt(data.num);