

%slfmriPlotTuningCurve.m
%
%
% author: steeve laquitaine
%purpose: plot tuning curve responses for different values of x0


function slfmriPlotTuningCurve(y,x0,colors)

st = SLmakeStat(y,x0);
var0Val = st.conditions(:,1);
myerrorbar(var0Val,st.mean,'yError',st.sem,...
    'Symbol=o',['Color=' num2str(colors)]);
%get legend
h(1) = plot(var0Val,st.mean,'.','linestyle','none','color','k');
xlabel(var0)
ylabel('Response')