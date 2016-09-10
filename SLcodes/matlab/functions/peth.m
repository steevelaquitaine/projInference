function [ output_args ] = peth( event2,bin,Xmin2,Xmax2 )
% Steeve laquitaine 09/06/2011

% PETH Summary of this function goes here
% plot a perievent histogram

%   Detailed explanation goes here
% creates perievent histogram.
% sp is a structure file containing timestamps for each trials.
% bin: the unit dividing the x axis (0.02 s for firing rate)
% Xmin and Xmax (of the graph)
% e.g.,
%  bin=0.02;
%  sp=a;
%  Xmin=b;
%  Xmax=c;
sp=event2.data;
Xmintmp=Xmin2.data;
Xmaxtmp=Xmax2.data;
edges=[Xmintmp:bin:Xmaxtmp];
psth=zeros(1,size(edges,2));
trl_nb=size(sp,2);
for j=1:trl_nb
    psth=psth+histc(sp(j).times,edges);
end
hold all
bar (edges,psth)
plot([0 0],[0 trl_nb],'r-');
xlim([Xmintmp Xmaxtmp])
ylim([0 max(psth)]);
xlabel({'Time unit(sec)' num2str(bin)});
ylabel('Firing rate (sp/s)');
text(-0.5,-max(psth)/6.5,event2.name,'fontweight','bold') % event
text(Xmax2.data-0.50,-max(psth)/6.5,Xmax2.name,'fontweight','bold') % axis' max
text(Xmin2.data-0.50,-max(psth)/6.5,Xmin2.name,'fontweight','bold') % axis' min
end

