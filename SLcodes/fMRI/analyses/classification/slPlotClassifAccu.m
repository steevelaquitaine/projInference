

%slPlotClassifAccu.m

%[160105]
%Classification
%Fisher, leaveOneOut

close all
figure('position',[231 425 958 319])

%Coherence
acc = [77.7 74.8 70.6 76.2 67.5 70.9 69.5 69.5 69.8 67.5 74.8 69.5 63 66.7 63.9 53.6 61.4 64];
ste = [2 2 2 2 2 2 2 2 2 2 2 2 3 2 3 3 3 3];
subplot(1,3,1)
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(50,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','rLO1','rLO2','IT','V7','IPS','BA39','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode coherence (cv ste)')


%Direction
acc = [29.7 27.3 26.3 24.6 21 26.6 28.6 27 27.7 31.6 29 28.6 25 26 26 21 25.6 26];
ste = [3 3 3 2 2 3 3 3 3 3 3 3 3 3 3 3 2 3 3];
subplot(1,3,2)
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(20,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','rLO1','rLO2','IT','V7','IPS','BA39','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode directions (cv ste)')


%Switching
acc = [67.9 65 57 64 54.7 57 48 59.7 55.8 60.9 60.9 58.9 53.1 56.6 55.8 53.9 53.5 52.7];
ste = [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];
subplot(1,3,3)
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(50,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','rLO1','rLO2','IT','V7','IPS','BA39','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode switching (cv ste)')


%% Switching removing coherence confound (by coh)
%because subjects switch more to prior at low coherence, trials tagged as
%switch to prior are more often low coherence trials and we might be
%decoding coherence instead of switch-to-prior
%load data with "slloadClassifAccu.m"
figure('color','w')
acc = [65 62 57 59 50 57 56.7 57 52 57.7 56.7 52.9 53.8 52.9];
ste = [3   3  3  3  3  3    3  3  3    3    3   3     3    3];

%coh0.06
subplot(1,2,1)
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(50,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode switching (cv ste) - coh 0.06')

%coh0.12
acc = [70 62 58 70   58   60   62   70  66   64   40   34   46 46];
ste = [6  6 6.9  6  6.9  6.9  6.8  6.4  6.6  6.7  6.9  6.6  7   7];
subplot(1,2,2)
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(50,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode switching (cv ste) - coh 0.12')


%% Switching removing coherence confound (by coh)
%10 - fold cross-validation
%load data with "slloadClassifAccu.m"
figure('color','w')
acc = [64 63 57 59 54 59 57 61];
ste = [3 3 4 3 4 3 4 3];

%coh0.06
subplot(1,2,1)
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(50,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode switching (10-fold cv, ste) - coh 0.06')

%coh0.12
acc = [60 50 38 55 55 57 63 63];
ste = [8 8 8 8 8 8 8 8];
subplot(1,2,2)
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(50,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode switching (10-fold cv, ste) - coh 0.12')


%% decoding directions when switch-to-prior
%Subjects might switch to prior because evidence is flat (Bayes) in which
%case we should not be able to decode motion directions
%note: theoretical chance accuracy at 1/4 = 0.25 not .2 because we don't
%use trials with 225 deg motion direction (because switching cannot be defined)
%#number instances for each class
%        sessions  total
%15  --> 3  2  3    8
%85  --> 5  2  2    9
%155 --> 22 8  15  45
%225 --> 0  0  0    0
%295 --> 22 5  15  52
figure('color','w')

%coh0.06
acc = [28    43  21.8   25   18 31  31  28  25  28  34  12  31 25];
ste = [7.9  8.7     7   7.6  7   8   8   8  7    8   8   6   8  7];
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(25,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode directions when subjects chose prior (cv ste) - coh 0.06')

%coh 0.12
%Cancelled analysis: too few samples
%#number instances for each class
%        sessions  total
%15  --> 1  0  1    2
%85  --> 0  0  0    0
%155 --> 6  0  4   10
%225 --> 0  0  0    0
%295 --> 8  0  5   13

%% decoding directions when switch-to-directions
%more of a control to check that oss of accuracy when subjects switched to
%prior is not due to smaller sample size but is a real effect. Here we
%should get good accuracy
%        sessions  total
%15  --> 11  10  13   34
%85  --> 11  10  13   34
%155 --> 17  12  20   49
%225 -->  0   0   0    0
%295 --> 17  17  21   55
figure('color','w')

%coh0.06
acc = [35 32 37.5  30 28.9 29.6 32 34 14.8 34 22.6  20  30  27];
ste = [ 4  4    4  4    4    4   4  4  3   4   3.6  3.5  4  3.9];
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(25,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode directions when subjects chose direction (cv ste) - coh 0.06')

%% decoding directions when switch-to-prior (pooled coh)
%Coherence and directions have no reason to be correlated when chosing
%prior so we pool all coherences to get more samples
%It doesn't help much. It just produce 9 instance/class 1 more than with 6%
%coherence.
%#number instances for each class
%        sessions  total
%15  -->  4   2  4  10
%85  -->  5   2  2   9
%155 -->  28  8  19  55
%225 -->  0   0  0    0
%295 -->  30  5  20  55
figure('color','w')

%coh0.06
acc = [28,33,19,36,33,25,36,22,42,33,28,25,33,22];
ste = [7,8,7,8,8,7,8,7,8,8,7,7,8,7];
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(25,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode directions when subjects chose prior (cv ste) - pooled coh')


%% decoding two directions (15 and 85) when switch-to-prior (low coherence (6%))
figure('color','w')
acc = [56];
ste = [12];
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(50,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'outsideBrain02','V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode directions when subjects chose prior (cv ste) - 6% coh')

%% decoding two directions (15 and 85) when switch-to-direction (low coherence (6%))
figure('color','w')
acc = [42];
ste = [6];
bar(acc,'facecolor',[.8 .8 .8],'edgecolor',[.5 .5 .5])
myerrorbar(1:length(acc),acc,'yError',ste,'Symbol=.','MarkerSize=1')
hline(50,'--r')
set(gca,'xtick',1:length(acc),'xticklabel',{'outsideBrain02','V1','V2','V3','V3A','MT','hV4','V7','IPS','FEF','dlPFC','vlPFC','vmPFC','OFC','aPFC'})
rotateXLabels(gca(),60)
axis tight
ylim([0 100])
xlim([0 length(acc)+1])
box off
title('Decode directions when subjects chose prior (cv ste) - 6% coh')




