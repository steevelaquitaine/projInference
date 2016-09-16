
%errorarea.m
%
%
% steeve laquitaine 12/03/2010 start: 17:30 ; end :18:17
% Debug1 16.03.2010; end 11:44
%
%
%usage:
%
%       y = errorarea(A,meancolor,fillcolor)
%
% INPUT: a matrix (average the column)
% OUTPUT: give an error area around the mean line

function y = errorarea(A,meancolor,fillcolor)

%checking
if all(isnan(A(:)))
    fprintf('%s \n','(errorarea) Error : The input matrix is all NaN');
    dbstack
    keyboard
end

moyenne=nanmean(A,2);
inan=find(isnan(moyenne)==1);
moyenne(inan)=[];
x=1:size(moyenne,1);
x=x';
% sem
for i=1:size(A,1)
    error(i,1) = std(A(i,:))/sqrt(length(A(i,:))-1);
end    
error(inan)=[];
errorupper=moyenne+error;
Eup=errorupper;
errorlower=moyenne-error;
Edn=errorlower;
%    meancolor=[0 1 0];
%    fillcolor=[0.8 1 0.8];
fill([x;x(end:-1:1)],[Eup;Edn(end:-1:1)],fillcolor,'edgecolor',[1 1 1])
hold on; plot(x,moyenne,'color', meancolor,'linewidth',1)