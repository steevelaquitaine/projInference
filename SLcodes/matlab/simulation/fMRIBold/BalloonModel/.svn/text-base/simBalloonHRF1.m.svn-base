% This is to simulate using balloon model

% Refer to equation 11 in 
% Mildner et al (2001) A qualitative test of the balloon model for BOLD-based MR signal changes at 3T. Magnetic resonance in medicine 46 (5) 891-9

% parameters
E0 = 0.4;
RT = 5;
DUR = 5;
dCBF = 0.7;
stimTime = [10];%ceil(rand(1,40)*100);
% stimTime = [
%    42.2000
%    76.1000
%   114.8000
%   150.0000
%   186.5000
%   222.9000
%   259.0000
%   296.3000
%   332.3000
%   365.6000
%   399.9000
%   437.1000
%   471.1000
%   506.9000
%   541.6000
%   575.7000
%   610.2000
%   645.9000
%   682.9000
%   719.2000];
stimTime = round(stimTime);
%stimTime = [100];
tao0 = 2;
taov = 30;
alpha = 0.4;

T = 757.2;

% fin kernel
fink = [[0:RT]/RT ones(1,DUR) [RT:-1:0]/RT];
fink = fink*(dCBF);

% simulate fin
%fin = ones(1,T);
%fin(stimTime:stimTime+length(fink)-1) = fin(stimTime:stimTime+length(fink)-1) + fink;
fin = zeros(1,T);
fin(stimTime) = 1;%rand(size(stimTime));
fin = conv(fin, fink) + 1;

%figure;plot(fin);

% now simulation

Mu = struct;
Mu.fin = fin;
Mu.E0 = E0;
Mu.tao0 = tao0;
Mu.taov = taov;
Mu.alpha = alpha;

tspan = [1,T];
y0 = [1; 1; 1];
ode = @(t,y) balloon_ode(t,y,Mu);
[t,y] = ode45(ode, tspan, y0);

% q = ones(1,T);
% v = ones(1,T);
% 
% for t=1:T
%     E(t) = 1-(1-E0).^(1/fin(t));
%     dq = fin(t)/tao0 * (E(t)/E0 - q(t)/v(t)) + 1/taov * (fin(t) - v(t).^(1/alpha))*q(t)/v(t);
%     dv =  1/taov * (fin(t) - v(t).^(1/alpha));
%     q(t+1) = q(t) + dq;
%     v(t+1) = v(t) + dv;
% end

%figure;plot(q); hold on;plot(v-q,'r')
%figure;plot(t,y(:,1)); hold on;plot(t, y(:,2)-y(:,1),'r')

tt = [1.1:0.1:T];
hbr = interp1(t,y(:,1),tt) - 1;
%hbo = interp1(t,y(:,2)-y(:,1),tt);
hbo = interp1(t,y(:,3)-y(:,1),tt);
v = interp1(t,y(:,2),tt);

% add noise?
noise = randn(size(hbr))/30 * 0;
noise2= randn(size(hbr))/60 * 0;
hbr = hbr + noise + noise2;
hbo = hbo - noise + noise2;
% filter
%hbr = passfilter(hbr, [0.5 0.01], 10);
%hbo = passfilter(hbo, [0.5 0.01], 10);
%%
figure;

plotTraces([hbo'-0.5 hbr']*2, [1 2], stimTime*10, 'rbk');
c = runningCorrelation(hbr,hbo, 40);
hold on;
plot(c,'k')
xlim([1 300])


%% all plots
figure;
plot(tt,hbr,'b');
hold on;
plot(tt,hbo,'r');
hold on;
plot(tt,v+1,'k');
hold on
fin = interp1([1:length(fin)], fin, tt);
plot(tt, fin+3, 'color', [0 0.5 0]);
hold on

fout = 1/(1+tao0/taov)*(tao0/taov*v.^(1/alpha) + fin);
plot(tt, fout+2, 'g');
hold on
e = 1-(1-E0).^(1./fin);
plot(tt, e+1, 'color', [.7 .7 .7]);

legend({'HbR','HbO','volume','flow in','flow out', 'e'})
%xlim([80 220])
xlim([0 80])
set(gcf,'color','w')
axis off