
% vmPdfs.m
%
%     author: Steeve Laquitaine
%       date: 150916
%    purpose: create von Mises circular distributions
%
%
%     inputs: 'norm' scales to probabilities, otherwise [];
%           when k are same, u must be a value of x for the code to work.
%
%
%      usage:
%
%           mPdfs = vmPdfs(1:1:360,[225 100],[0 10],'norm')
%
%
%Description:
%
%       vm = exp(k*cos(w*(xrad-urad(1)))-k) / (2*pi.*besseli(0,k,1))
%
%   This code produces von Mises densities vm(u,k) based on the
%   equation vm=exp(k.*cos(x-u))./(2*pi.*besseli(0,k)); The code works for any
%   value of k (but not for inf).The equation is adjusted because of the
%   following numerical issues: when k>700, vm is NaN because besseli(0,k) and
%   exp(k.*cos(x-u)) reach numerical limits. exp(k.*cos(x-u)-k) scales vm
%   without changing its shape. besseli(0,k,1)) does same. The adjusted
%   equation and the exact equation yield exact same results except that the
%   adjusted equation works for large k (>>700).
%
%   When k > 1e300 densities tend toward a delta density. Due to numerical
%   limits matlab outputs NaN so we generate delta density for k>300.

function mPdfs = vmPdfs(x,u,k,type)


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

%radians
xrad = SLde2r(x,1); xrad=xrad';
urad = SLde2r(u,1); urad=urad';
k=k';
%When von mises with different mean u1,u2,u3 but with same k are input
%We can get von Mises with mean u2,u3,etc...simply by rotating the von
%mises with mean u1 by u2-u1, u3-u1 etc...
%When we don't do that we get slightly different von mises with different
%peakvalue due to numerical instability caused by cosine and exponential
%functions.
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
        
        %-------------------------
        %case k tends toward + inf
        %-------------------------
        %delta function
        if k(1) > 1e300
            
            mPdfs = zeros(length(xrad),1);
            mPdfs(xrad == urad(1)) = 1;
            %make other densities by circshifting the first
            for i = 2 : numel(u)
                rotation = find(x==u(i)) - find(x==u(1));
                mPdfs(:,i) = circshift(mPdfs(:,1),[rotation,0]);
            end
            
        else
            
            %---------
            %k non inf
            %---------
            w = 1;
            k = k(1);
            
            %should work with parallel computing (parfor loop). ".*" seems to
            %work badly with parfor. But I think k must be unique.
            mPdfs = exp(k*cos(w*(xrad-urad(1)))-k)./(2*pi.*besseli(0,k,1));
            %make other densities by circshifting the first
            for i = 2 : numel(u)
                rotation = find(x==u(i))-find(x==u(1));
                mPdfs(:,i) = circshift(mPdfs(:,1),[rotation,0]);
            end
            
            %scale to probabilities.
            if strcmp(type,'norm')==1
                %mPdfs=mPdfs./sum(mPdfs(:,1));
                mPdfs=mPdfs/sum(mPdfs(:,1));
            end
        end
    end
end

%When von mises with different mean and k are input, calculate densities
%separately. Make sure each mean u maps with a k.
if sum(k - k(1))~=0
    
    %case k not too high
    k     = k(ones(numel(xrad),1),:);
    x2rad = xrad(:,ones(numel(urad),1));
    u2rad = urad(ones(numel(xrad),1),:);
    k2    = k(ones(numel(xrad),1),:);
    w     = 1;
    mPdfs = exp(k2.*cos(w*(x2rad-u2rad))-k2)./(2*pi.*besseli(0,k2,1));
    
    %case k tends toward inf, delta density
    kinfCols      = find(k2(1,:)>300);
    u2radkinfCols = u2rad(:,kinfCols);
    deltas        = zeros(numel(xrad),length(kinfCols));
    deltas(x2rad(:,kinfCols)-u2radkinfCols==0) = 1;
    mPdfs(:,kinfCols) = deltas;
    
    %for parallel processing (bu 30 times slower)
    %     for i=1:size(k2,1)
    %         for j=1:size(k2,2)
    %             mPdfs(i,j)=exp(k2(i,j)*cos(w*(x2rad(i,j)-u2rad(i,j)))-k2(i,j))/(2*pi*besseli(0,k2(i,j),1));
    %         end
    %     end
    
    %scale to pdfs.
    if strcmp(type,'norm')==1
        Z_=sum(mPdfs);
        Z=Z_(ones(numel(x),1),:);
        mPdfs=mPdfs./Z;
    else
    end
end
%degrees to radians

function radians = SLde2r(ang,sign)

%call for help (may slow down code)
% if ieNotDefined('ang')
%     help de2r
%     return
% end

%[150907]
% %not signed radians (1:2*pi)
% radians = (ang/360)*2*pi;
%
% %sign radians(-pi:pi)
% if sign==1
%     radians(ang>180)=(ang(ang>180)-360)*(2*pi/360);
% end

%memory and speed
ang   = single(ang);
twopi = single(2*pi);

%not signed radians (1:2*pi)
radians = (ang/360)*twopi;

%sign radians(-pi:pi)
if sign==1
    radians(ang>180)=(ang(ang>180)-360)*(twopi/360);
end

radians = double(radians);
