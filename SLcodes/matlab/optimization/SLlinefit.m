
%SLlinefit.m
%
% author: Steeve Laquitaine
%   date: 140820
%purpose: linear fit to data
%
% usage:
%
%   SLlinefit(1:10,[1:10]+randn(1,10),'r');
%
%varargin
%--------
%   'noNegativeValue': show only predictions > 0.
%
%
%
%examples
%---------
%('noNegativeValue')
%
%  SLlinefit(1:10,[1:10]+randn(1,10),'r','noNegativeValue');
 
function [xfit,yfit] = SLlinefit(x,y,linecolor,varargin)

%available data
i =~ isnan(y);

%case there is no data
if sum(i)==0
    xfit = NaN(length(y),1);
    yfit = NaN(length(y),1);
    return
end

y = y(i);

%case default
%------------
if isempty(varargin)

    %fit
    xfit = x;
    if size(x,1)~=size(y,1)
        x = x';
    end
    
    %fit
    P = polyfit(x(i),y,1);
    yfit = polyval(P,x);
    
    %plot
    plot(xfit,yfit,'-',...
        'linewidth',0.5,...
        'color',linecolor,...
        'LineSmoothing','on');
end


%case we do not want negative values
%------------------------------------
if strcmp(varargin,'noNegativeValue')
    
    xi = x(i);
    xi_contSpace = xi(1):0.01:xi(end);
    
    %fit linear weights P
    P = polyfit(xi,y,1);
    
    %get best linear predictions
    xfit = xi_contSpace;
    yfit = polyval(P,xi_contSpace);
    
    %plot positive values
    plot(xi_contSpace(yfit>0),yfit((yfit>0)),'-',...
        'linewidth',0.5,...
        'color',linecolor,...
        'LineSmoothing','on');    
end