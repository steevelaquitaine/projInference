
 %author: steeve laquitaine
   %date: 140310
%purpose: z-score data: normalization in std to the mean
  %usage: x=SLzscore(x,'whole')

function x=SLzscore(x,type)
%mean and std are computed over whole matrix
if strcmp(type,'whole')
    x=(x-mean(x(:)))/std(x(:));
end
%mean and std are computed over columns
if strcmp(type,'col')
    MeanX=mean(x,1);
    StdX=std(x);
    x=(x-MeanX(ones(size(x,1),1),:))./StdX(ones(size(x,1),1),:);
end
%mean and std are computed over raws
if strcmp(type,'raw')
    MeanX=mean(x,2);
    StdX=std(x')';
    x=(x-MeanX(:,ones(size(x,2),1),:))./StdX(:,ones(size(x,2),1),:);
end