function [tp_to_evt, mintp_to_evt,mxtp_to_evt] = raster(unit_times, event, Xmin, Xmax)
%PETH Summary of this function goes here
%   Detailed explanation goes here

% % e.g
% unit_times=unit_01_flag1;
% ITI=2; %  2-2.5 sec
% Xmin=S_trial_begin;
% Xmax=Reinforcement+ITI;
% event=Reinforcement;

Xmintmp=Xmin.data;
Xmaxtmp=Xmax.data;

% sort spikes in trials
trl_nb=size(event.data,1);
for i=1: trl_nb
    [~,spike(i).times]=find_betw(unit_times, Xmintmp(i),Xmaxtmp(i),'strict');
end

% raster plots
clear i
hold all;
tt0=event.data;
tp_to_evttp=[];
for t=1:trl_nb % trial (y-axis)
    tt=spike(t).times;
    hgt(t)=t;
    for i=1: length(tt) % time (x-axis)
        tp_to_event(i)=tt(i)- tt0(t);
        line([tp_to_event(i)  tp_to_event(i)], [hgt(t)-1 hgt(t)]);
    end    
    tp_to_evttp=[tp_to_evttp; tp_to_event'];
    tp_to_evt(t).times=tp_to_event;
    clear tt
    clear tp_to_event
end
% axis
mintp_to_evt=min(min(tp_to_evttp));
mxtp_to_evt=max(max(tp_to_evttp));

xlabel('time(s)');
ylabel('trials');
ylim([0 trl_nb]);
xlim([mintp_to_evt mxtp_to_evt]);
plot([0 0],[0 trl_nb],'r-');
text(-0.5,-50,event.name,'fontweight','bold') % event
text(mxtp_to_evt-0.50,-trl_nb/6.5,Xmax.name,'fontweight','bold') % axis' max
text(mintp_to_evt-0.50,-trl_nb/6.5,Xmin.name,'fontweight','bold') % axis' min
end

