
close all

x=smooth_nolearning
nuance=10

% Enter parameter
n =size(x,2)+nuance % must be >= to the number of ploted curve


for i=1:size(x,2)
    map = colormap(gray(n));
    hold on;
    plot1=plot(1:size(x,1),x(:,i),'-','color',map(i,:),'linewidth',3);
end
