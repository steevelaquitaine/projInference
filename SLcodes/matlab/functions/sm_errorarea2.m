
function [ y ] = sm_errorarea(A,meancolor,fillcolor);


% steeve laquitaine 12/03/2010 start: 17:30 ; end :18:17
% Debug1 16.03.2010; end 11:44

% INPUT: a matrix (average the column)
% OUTPUT: give an error area around the mean line


wind=11; % frame length
degree=1;

moyenne=nanmean(A,2);
inan=find(isnan(moyenne)==1);
moyenne(inan)=nanmean(moyenne);
x=1:size(moyenne,1);
x=x';

% sem
for i=1:size(A,1)
    error(i,1)=sem(A(i,:));
end    
inan=find(isnan(error)==1);
error(inan)=nanmean(error);

errorupper=moyenne+error;
Eup=sgolayfilt(errorupper,degree,wind);
errorlower=moyenne-error;
Edn=sgolayfilt(errorlower,degree,wind);

%    meancolor=[0 1 0];
%    fillcolor=[0.8 1 0.8];

fill([x;x(end:-1:1)],[Eup;Edn(end:-1:1)],fillcolor,'edgecolor',[1 1 1])
hold on; plot(x,sgolayfilt(moyenne,degree,wind),'color', meancolor,'linewidth',1)