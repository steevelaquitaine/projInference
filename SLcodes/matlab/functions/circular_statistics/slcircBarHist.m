


%slcircBarHist(1:1:360,10)


function slcircBarHist(x,bin)

keyboard

binStart = 1:bin:360;
binEnd = bin:bin:360;

for i = 1 : length(binStart)    
     binCenter = binStart(i) + bin/2;
     %find values within bin distance to bin center
     maxDis = SLvectors2signedAngle(binCenter,binEnd(i),'polar');
     minDis = SLvectors2signedAngle(binCenter,binStart(i),'polar');     
     xtoCenter = SLvectors2signedAngle(binCenter,x,'polar');
     
     find(xtoCenter > )
     
     
%     count(i) = length(find(d20 > dbinStart(i) & d20 <= dbinEnd(i)));
end

bar(1:bin:360,count)

