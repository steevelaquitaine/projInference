
%slmakeGaussianPdfs.m
%
%      author: steeve laquitaine
%        date: 141105
%     purpose: compute a 1 dimensional (x) Gaussian
%
%       usage:
%
%              p = slmakeGaussianPdfs(-100:1:100,[0 1],[1 2],'norm')
%
%
%Description:
%
%       x = space
%       u = mean of x
%       k = stdev of x
%       type: 'norm' if you want probability that
%              sum to 1.
%
%
function pdf = slmakeGaussianPdfs(x,u,k,type)

%check that x is a row vector
if size(x,1)>size(x,2)
    x=x';
end

%check that u and k are col vectors
if size(u,1)<size(u,2)
    u=u';
end
%check that u is a col vector
if size(k,1)<size(k,2)
    k=k';
end

%case all k are same
if sum(k - k(1))==0
    %if mean u is not one of x
    if isempty(intersect(x,u)==0)
        fprintf('\n %s \n','(vmPdfs) WARNING : Mean has to be one of x values. Change mean.....')
        if isnan(u)
            sprintf('(vmPdfs) WARNING : Von Mises Mean is NaN....')
        end
        dbstack
        keyboard
    else
        
        %Gaussian (amplitude=1)
        for i = 1 : length(u)
            dist = x - u(i);
            pdf(:,i) = exp(-0.5*(dist./k).^2);
        end
        
        %scale to probabilities.
        if strcmp(type,'norm')==1
            pdf = pdf*1./(sqrt(2*pi)*k);
        end
    end
end

%case k are different
if sum(k - k(1))~=0
    
    %case k not too high
    k     = repmat(k',length(x),1);
    x     = repmat(x',1,length(u));
    u     = repmat(u',size(x,1),1);
    dist  = x - u;
    pdf   = exp(-0.5*(dist./k).^2);    
    
    if strcmp(type,'norm')==1
        %scale to probabilities.
        if strcmp(type,'norm')==1
            pdf = pdf*1./(sqrt(2*pi)*k);
        end
    end
end

