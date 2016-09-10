
%d
%d are directions in radians (corresponding degrees range from 0:360).

%k:  
%k is the precision or concentration parameter (1/width). It is 42 in
%Girshick et al,2011, for orientations because it produces a biologically 
%plausible population average tuning width of 17 degrees. 
%Girshick, A. R., Landy, M. S. & Simoncelli, E. P. Cardinal rules:
%visual orientation perception reflects knowledge of environmental 
%statistics. Nature Publishing Group 14, 926?932 (2011).
%Swindale, N. V. Orientation tuning curves: empirical description 
%and estimation of parameters. Biol Cybern 78, 45?56 (1998).


%ui
%ui is the preferred orientation (the phase shift of the cosine).
%Swindale, N. V. Orientation tuning curves: empirical description 
%and estimation of parameters. Biol Cybern 78, 45?56 (1998).

%u=3:6:360 to be sure one of the tuning curve is centered on the prior
%mean.


%w
%w is the frequency. It is 2 in Girshick's paper because orientation
%selectivity peak every 180 degrees that is two times in 360
%degrees. In the case of motion direction w=1.

%The mean of a Von Mises is its mode(easier to calculate) but if the mean
%is not a contained in the sample d (e.g., 0.5 while d is 1:1:360) the mode
%and is an approximation.

%%compute weighted sum of cos and sin of angles 
%r = w'*exp(1i*alpha);
%
%%obtain mean by 
%mu = angle(r);
%http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

%alpha:angle in radian


%References
%http://mathworld.wolfram.com/vonMisesDistribution.html
%Girshick, A. R., Landy, M. S. & Simoncelli, E. P. Cardinal rules:
%visual orientation perception reflects knowledge of environmental 
%statistics. Nature Publishing Group 14, 926?932 (2011).
%Swindale, N. V. Orientation tuning curves: empirical description 
%and estimation of parameters. Biol Cybern 78, 45?56 (1998).
%http://www.mathworks.com/matlabcentral/newsreader/view_thread/237503

%Population tuning curve
function TC=tuningC
%degrees (easier)
deg=1:1:360;
u=3:6:360;
k=42; 
for i=1:60
    TC.i(i,:)=tuningCi(deg,u(i),k);
end

%population average(not sure about the operation)
TC.av=mean(TC.i);
hold on; plot(deg,TC.av,'k','linewidth',2)

%A neuron's tuning curve
function TC=tuningCi(deg,ui,k)
w=1;

%degrees to radian
d=de2r(deg);
ui=de2r(ui);

%tuning curve
%Z=sum(exp(k*cos(1*(d-ui))));
TC=exp(k*cos(w*(d-ui)))/2*pi*besseli(0,k);
plot(deg,TC,'-','color',[.3 .3 .3]);

%check mean and std (wrong!!!)
%TC.u=ra2d(sum(TC.i.*d));
%TC.s=sqrt(1-(besseli(1,k)/besseli(0,k)));
%TC.s2=sqrt(sum(TC.i.*(d-225).^2)/sum(TC.i));

%graphics
hold on;
xlim(de2r([min(deg) max(deg)]))
xlim([min(deg) max(deg)])

ylim([0 max(TC)])
% set(gca,'xtick',d(1:10:numel(deg)),'xticklabel',deg(1:10:numel(deg)))
ylabel('Firing rate (arbitrary unit)')
xlabel('Motion direction (degree)')
box off





%Nested
%convert degrees to radians
function radians = de2r(angle)
radians = (angle/360)*2*pi;

%convert radians to degrees
function degrees = ra2d(angle)
degrees = (angle/(2*pi))*360;

% if larger than 360 degrees then subtract
% 360 degrees
while (sum(degrees>360))
  degrees = degrees - (degrees>360)*360;
end

% if less than 360 degreees then add 
% 360 degrees
while (sum(degrees<-360))
  degrees = degrees + (degrees<-360)*360;
end
