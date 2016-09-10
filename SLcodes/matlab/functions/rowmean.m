function y=rowmean(x)

%steeve laquitaine 12042009

%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

%  use a matrix
%  output :Standard error of the mean std/sqrt(count) for each line

for i=1:length(x(:,1))
    y(i,1)=nanmean(x(i,:))
end
    