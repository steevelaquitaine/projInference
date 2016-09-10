
%Usage
    %k.m=42; 
    %k.p=42; 
    %[mPdfs,l,pr,po,di,m,MAPfull,MAPpdfs]=GirshickEstimator(k);

%di
%d are directions converted in radians (corresponding degrees range from 
% 0:360). di should always be 1:1:360 because I saw patterns I can not
% explain if di is for example 5:1:355 at llh mean=45 degrees.

%k:  
%k is the precision or concentration parameter (1/width). It is 42 in
%Girshick et al,2011. k is a fit parameter.
%*IMPORTANT* "In our case measurement distributions and likelihoods are identical because
%k is assumed to be same for all motion directions."
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

%Likelihood
%The likelihood is the transpose of the measurement
%distribution that is then normalized.

%Girshick Bayesian estimator
%Each single trial's measurement is transformed into a map estimate 
%(argmax posterior). This mapping between measurement and estimates is the
%estimator.


%Test the behavior of the model
%................................
%Notes about the model:
%-When prior width = measurement width, MAPs are evenly distant by
%0.5 degrees. MAPs are independent on the size of the width and
%independent on the resolution of directions.

%MAPs distance only depends on measurement resolution. 
%It increases as measurement resolution decreases.

%for a given range of widths [1 40] or [40 1], I can get MAPs that all 
%differ but not in smooth steps so that I can match any given data to a
%MAP. Resolution of m=2; resolution of di=0.0001.


%Test the behavior of the model
%................................
%k.m=1; k.p=40 is the case where posterior virtually equals the prior.
%k.m=40;k.p=1  is the case where posterior virtually equals the prior.
%Are there unique pairs k.m and k.p that explains the data in each 
%conditions ? I think so. I think each pair predicts a unique combination 
%of variability and mean of the data. 


%References
%http://mathworld.wolfram.com/vonMisesDistribution.html
%Girshick, A. R., Landy, M. S. & Simoncelli, E. P. Cardinal rules:
%visual orientation perception reflects knowledge of environmental 
%statistics. Nature Publishing Group 14, 926?932 (2011).
%Swindale, N. V. Orientation tuning curves: empirical description 
%and estimation of parameters. Biol Cybern 78, 45?56 (1998).
%http://www.mathworks.com/matlabcentral/newsreader/view_thread/237503

function [mPdfs,l,pr,po,di,m,MAPfull,MAPpdfs,means,stds]=GirshickEstimator(k)
%direction(degrees)
di=[1:1:360]';
%di=fix(di*10)/10;%2 floating points

%measurements(degree)
m=[1:1:360];
n=numel(m);

%measurement distributions(pdf)
mPdfs=vmPdfs(m,di,k.m);

%likelihood(rows)
l=mPdfs;
Zl=sum(l,2);
Zl=Zl(:,ones(numel(di),1));
l=l./Zl;

%priors(rows)
pr.m=ones(n,1)*225;
pr.pdf=vmPdfs(di',pr.m,k.p);
pr.pdf=pr.pdf';

%posteriors(rows)
po=l.*pr.pdf;
%po=bsxfun(@rdivide,po,sum(po,2));
Zpo=sum(po,2);
Zpo=Zpo(:,ones(numel(di),1));
po=po./Zpo;
po=fix(po*10^10)/10^10;%avoid computing errors at <10^10 floating points;

%MAP estimates
[X,I]=max(po');
MAP=di(I);

%MAP when measurement=45
clear I
[X,I]=max(l(m==45,:)');
MAP(m==45)=di(I);


%set probability of MAPs never observed given a direction to zeros. This
%doesn't work need to do sthg else.
%MAP0=[];%setdiff(di,MAP);
%MAPfull=[MAP;MAP0];
%MAPpdfs=[mPdfs;zeros(numel(MAP0),numel(di))];
%for now this is better
MAP0=[];%setdiff(di,MAP);
MAPfull=[MAP;MAP0];
MAPpdfs=[mPdfs;zeros(numel(MAP0),numel(di))];



%Plots
%plot matrices
%drawGirshickInference(di,m,l,pr,po)
%hold on; 
%plot(MAP,m,'.r','Markersize',6); 
%xlim([0 360]); 
%ylim([0 360])
%plot slices of the matrices
%drawSlices(di,m,l,pr,po)
%plot Pdfs of measurements and estimates given direction i.
%variability of MAP = variability of m, but the prior distorts MAP's pdf 
%compared to m's pdf.
[means,stds]=drawEstimatePdfs(di,m,mPdfs,MAPfull,MAPpdfs);





%Nested
%BLS estimation
function Av=vmBLSestimator(di,p)
%for test with a von Mises.
%u=180;
%k=42;
%p=vmPdf(di,u,k);


%MAP estimate
%convert directions into vectors
v=polar2cartesian(di,p);

%vectors are averaged.
Av=vectorAverage(v);
%Plots
function drawGirshickInference(di,m,l,pr,po)
figure('color','w'); 

subplot(1,3,1)
imagesc(di,m,l);
title('Measurement(colum slices) & Likelihood(row slices) distributions')
ylabel('Measurements(degrees)')
xlabel('Directions(degrees)')
set(gca,'YDir','normal');
colormap(gray);
axis square

subplot(1,3,2)
imagesc(di,m,pr.pdf);
title('Prior distribution')
xlabel('Directions(degrees)')
ylabel('Measurements(degrees)')
set(gca,'YDir','normal');
colormap(gray);
axis square

subplot(1,3,3)
imagesc(di,m,po);
title('Posterior distribution')
xlabel('Directions(degrees)')
ylabel('Measurements(degrees)')
set(gca,'YDir','normal');
colormap(gray);
axis square
function drawSlices(di,m,l,pr,po)
figure('color','w'); 

%examples llh.
exples=[find(m==45) 
find(m==90) 
find(m==135) 
find(m==180)
find(m==225)
find(m==270)
find(m==360)];
for i=1:numel(exples)
    subplot(numel(exples),1,i)
    hold all
    %llh
    plot(di,l(exples(i),:),'color','k')
    %prior
    plot(di,pr.pdf(exples(i),:),'--','color',[.5 .5 .5])
    %posterior
    plot(di,po(exples(i),:),'b','linewidth',1)
    text(di(di==m(exples(i))),-0.015,num2str(di(di==m(exples(i)))));
    text(pr.m(exples(i)),-0.015,num2str(pr.m(exples(i))));
    set(gca,'xtick',[],'xticklabel',[])
    box off
end
function [means,stds]=drawEstimatePdfs(di,m,mPdfs,MAPfull,MAPpdfs)
figure('color','w'); 

%examples directions.
ex=[45 90 135 180 225 270 360];
exples=[find(di==ex(1)) 
find(di==ex(2)) 
find(di==ex(3)) 
find(di==ex(4))
find(di==ex(5))
find(di==ex(6))
find(di==ex(7))];
nx=numel(exples);

%circular analyses
coord=polar2cartesian(MAPfull,1);

for i=1:numel(exples)
    %measurement distribution
    subplot(3,nx,i);
    plot(m,mPdfs(:,exples(i)),'.','color',[.5 .5 .5],'linewidth',2);
    xlim([0 360])
    box off
    title(strcat(num2str(di(exples(i))),' degrees'))
    ylabel('Probability')
    xlabel('Measurement')
    
    %estimate distribution
    subplot(3,numel(exples),i+nx);
    hold all
    plot(MAPfull,MAPpdfs(:,exples(i)),'.r','linewidth',2)
    plot([225 225],[0 max(MAPpdfs(:,exples(i)))],'b--')
    text(225,max(MAPpdfs(:,exples(i)))+0.003,'prior')
    xlim([0 360])
    ylabel('Probability')
    xlabel('MAP estimate')
    box off
    
    %Circular mean & std of estimate distribution
    subplot(3,numel(exples),2*nx+i);
    data=CircWeightedMeanStd(coord,MAPpdfs(:,exples(i)));
    means=data.deg.mean;
    stds=data.deg.std;
    
    myerrorbar(ex(i),means,'yError',stds,'Color=[1 0 0]','Symbol=o')
    ylim([110 300])
    ylabel('Estimates mean & std (degrees)')
    xlabel('Displayed directions (degrees)')
end





%von Mises
function mPdfs=vmPdfs(x,u,k)
w=1;
%radian
u=de2r(u); u=u';
x=de2r(x); x=x';
u2=u(ones(numel(x),1),:);
x2=x(:,ones(numel(u),1));
%mPdfs=nan(numel(x),numel(u));
mPdfs=exp(k.*cos(w*(x2-u2)))/2*pi.*besseli(0,k);
%pdf
Z_=sum(mPdfs);
Z=Z_(ones(numel(x),1),:);
mPdfs=mPdfs./Z;
%vector average
function [deg rad cart]=vectorAverage(v)
%v is a directional vector (cartesian coordinates)
%r is the unit circle radius.

% %cartesian of the mean
% Av.cart=nanmean(v,1);
% x=Av.cart(:,1);
% y=Av.cart(:,2);
% 
% %radian (arctan(opposite/adjacent))
% %Av.rad=atan(y./x);
% [Av.rad, Av.deg]=cart2polar(Av.cart);

%cartesian of the mean
cart=nanmean(v,1);
x=cart(:,1);
y=cart(:,2);

%radian (arctan(opposite/adjacent))
%Av.rad=atan(y./x);
[rad, deg]=cart2polar(cart);
%polar to cartesian
function v=polar2cartesian(di,r)
%di:angs in degree
%r:radius of the unit circle

%degrees & radians
di2.deg=di;
di2.rad=(di2.deg*pi)/180;
% di2.rad=(di2.deg*3.14)/180;


%cartesian vectors
x=r.*cos(di2.rad);
y=r.*sin(di2.rad);
v=[x y];
%polar to cartesian
function [rad,deg]=cart2polar(v)
%v is the cartesian vector

%radian (arctan(opposite/adjacent))
x=v(:,1);
y=v(:,2);
rad=atan(y./x);

%degrees
deg=ra2d(rad);

%adjust each angle according to his quadrant (in degree) otherwise atan
%gives wrong values for negative x and y.
for i=1:numel(x)
    if x(i)>=0 && y(i)>=0              %(quadrant 1)
        deg(i)=deg(i);
    elseif x(i)<0                      %(quadrant 2 & 3)
        deg(i)=deg(i)+180;
    elseif x(i)>=0 && y(i)<0           %(quadrant 4)
        deg(i)=deg(i)+360;
    end
end
%degrees to radians
function radians=de2r(ang)
radians=(ang/360)*2*pi;
%convert radians to degrees
function degrees=ra2d(theta)
degrees=(theta/(2*pi))*360;

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
%polar to cartesian
function v=polar2cartesian(di,r)
%di:angs in degree
%r:radius of the unit circle

%degrees & radians
di2.deg=di;
di2.rad=(di2.deg*pi)/180;

%cartesian vectors
x=r.*cos(di2.rad);
y=r.*sin(di2.rad);
v=[x y];
%Calculate the circular statistics of the data
function data=statcircular(coord)
%vector of cartesian coordinates (coord)

%coordinates
data.coord.all=coord;

%from cartesian to angles (degree)
data.deg.all=getangle(coord(:,1),coord(:,2));

%mean's cartesian & degree
data.coord.mean=nanmean(coord,1);
data.deg.mean=getangle(data.coord.mean(:,1),data.coord.mean(:,2));

%calculate the std (degree)
%Apply the rule of thumb that follows. It seems to work fine intuitively. 
%initialize the 'sample' and 'mean' variables used to calculate the std
data.num=numel(data.deg.all);
data.deg.allforstd=data.deg.all;
data.deg.meanforstd=repmat(data.deg.mean,data.num,1);

%std
%mean is between 0 and 180.
if data.deg.mean+180<=360
    %if estimation is>=mean direction+180
    for i=1:data.num
        if data.deg.all(i)>=data.deg.mean+180
            data.deg.allforstd(i)=data.deg.all(i)-360;
        end
    end
    %mean is between 180 and 360.
else
    %if the estimated direction sampled is <=the mean direction - 180
    for i=1:data.num
        if data.deg.all(i)<=data.deg.mean-180
            data.deg.meanforstd(i)=data.deg.mean-360;
        end
    end
end

%var
data.deg.var=nanmean((data.deg.allforstd-data.deg.meanforstd).^2,1);

%std
data.deg.std=sqrt(data.deg.var);

%sem
data.deg.sem=data.deg.std/sqrt(data.num);
%Calculate the circular statistics of the data
function data=CircWeightedMeanStd(coord,p)
%vector of cartesian coordinates (coord)

%coordinates
data.coord.all=coord;

%cartesian to angles (degree)
data.deg.all=getangle(coord(:,1),coord(:,2));

%mean's cartesian & degree
data.coord.mean=sum(p(:,ones(1,2)).*data.coord.all);
data.deg.mean=getangle(data.coord.mean(:,1),data.coord.mean(:,2));

%std (degree).
data.num=numel(data.deg.all);
data.deg.allforstd=data.deg.all;
data.deg.meanforstd=repmat(data.deg.mean,data.num,1);

%mean is between 0 and 180.
if data.deg.mean+180<=360
    %if estimation is>=mean direction+180
    for i=1:data.num
        if data.deg.all(i)>=data.deg.mean+180
            data.deg.allforstd(i)=data.deg.all(i)-360;
        end
    end
    %mean is between 180 and 360.
else
    %if the estimated direction sampled is <=the mean direction - 180
    for i=1:data.num
        if data.deg.all(i)<=data.deg.mean-180
            data.deg.meanforstd(i)=data.deg.mean-360;
        end
    end
end

%var
data.deg.var=sum(p.*((data.deg.allforstd-data.deg.meanforstd).^2));

%std
data.deg.std=sqrt(data.deg.var);

%sem
data.deg.sem=data.deg.std/sqrt(data.num);
%Convert from cartesian coordinates to angles (degree)
function [angle]=getangle(x,y)
% check! to check if the function works fine, write
% e.g., input=180; output=getangle(cos(angle*pi/180),sin(angle*pi/180));
% if the function works input and output should always be the same between
% 0 and 360.

% convert from cartesian coordinates to angle in radians
angle = atan(y./x); % theta = arctan(opposite/adjacent);

% adjust each angle according to his quadrant (in degree)
for i=1:numel(x) % sample each angle
    if x(i)>=0 && y(i)>=0                   %(quadrant 1)
        angle(i) = angle(i)*180/pi;
    elseif x(i)<0                      %(quadrant 2 & 3)
        angle(i) = angle(i)*180/pi + 180;
    elseif x(i)>=0 && y(i)<0               %(quadrant 4)
        angle(i) = angle(i)*180/pi + 360;
    end
end


