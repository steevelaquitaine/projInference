

%SLdrawVonMisesStdofKarea.m
%
% author: steeve laquitaine
%   date: 140904
%purpose: plot a von Mises with with its average K parameter and its
%         error area calculated from std of K.
%  usage:
%
%       Explanation
%       %----------
%       SLdrawVonMisesStdofKarea(-180:1:180,0,3,4,'Explanation')
%
%
%       Analysis
%       %-------
%       SLdrawVonMisesStdofKarea(-180:1:180,0,3,4)

function [meanVonMises,ErrorKup,ErrorKdn] = SLdrawVonMisesStdofKarea(x,vmMean,...
    vmMeanOfK,stdOfK,varargin)

%mean von Mises
meanVonMises = vmPdfs(x,vmMean,vmMeanOfK,'norm');

%min and max values of K
Kup = vmMeanOfK + stdOfK;
Kdn = vmMeanOfK - stdOfK;

%K is defined >0
%K<0 has no meaning.
Kdn(Kdn<0)=0;

%Step 1. plot von Mises for example K values in between Kup and Kdn
%[Kdn:.1:Kup Kup] : repeat Kup to make sure it is counted.
RangeOfK = [Kdn:.1:Kup Kup];
errorDntoUp = vmPdfs(x,linspace(0,0,numel(RangeOfK)),RangeOfK,'norm');


%case explanation of the analysis
if sum(strcmp(varargin,'Explanation'))
    
    % %plot
    % set(gcf,'color','w')
    % subplot(1,2,1)
    % axis square
    % hold all
    % box off
    % plot(x,errorDntoUp,'color',[.7 .7 .7],'linesmoothing','on')
    % plot(x,meanVonMises,'color','k','linesmoothing','on')
    %
    % title('von Mises between k+std and k-std')
    %
    % %Step 2. plot von Mises error area between Kup and Kdn
    ErrorKdn = min(errorDntoUp');
    ErrorKup = max(errorDntoUp');
    % plot(x,ErrorKup,'color',[1 0 .0],'linesmoothing','on')
    % plot(x,ErrorKdn,'color',[0 0 1],'linesmoothing','on')
    %
    % subplot(1,2,2)
    % axis square
    % box off
    hold all
    SLerrorarea(meanVonMises,[],x,ErrorKup,ErrorKdn,'inputErrors');
    plot(x,ErrorKup,'--','color',[1 0 .0],'linesmoothing','on',...
        'linewidth',0.5)
    plot(x,ErrorKdn,'--','color',[0 0 1],'linesmoothing','on',...
        'linewidth',0.5)
    
    % title('Std area between k+std and k-std')
    
    %clean up
    % SLConventionUp(gcf)
else
    
    %case plot
    ErrorKdn = min(errorDntoUp');
    ErrorKup = max(errorDntoUp');
    hold all
    SLerrorarea(meanVonMises,[],x,ErrorKup,ErrorKdn,'inputErrors',...
        varargin{:});
end






