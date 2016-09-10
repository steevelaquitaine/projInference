
%   SLsimulateBAYESmapAndSamplingAndCompetitionPredictions.m
%
%     author: steeve laquitaine
%       date: 140601 last modification 140610
%
%    purpose: simulate motion direction estimation behavior with Bayesian
%    inference (learnt priors, cardinal priors, motor noise, random estimation).
%
%description: 
%
%   e.g., 1 fig --> 1 x 1 figures
%   e,g,, 2 fig --> 1 x 2 figures
%   e.g., 3 fig --> 1 x 3
%   e.g., 4 fig --> 2 x 2 
%   e.g., 5 fig --> 2 x 3 
%   e.g., 6 fig --> 2 x 3 
%   e.g., 7 fig --> 2 x 4 
%   e.g., 8 fig --> 3 x 3 
%   e.g., 9 fig --> 3 x 3 


function SLpositionFigures

%get the N figures opened
figH=SLgetFigures;

%divide the screen in N squares
numFig=numel(figH);

%Divide the screen in columns and rows
sqrtN=sqrt(numFig);
if SLisinteger(sqrtN)==1
    numCols=sqrtN;
    numRows=sqrtN;
else
    %otherwise
    numRows=fix(sqrtN);
    numCols=ceil(numFig/numRows);
end

%adjust by one figure
if numRows*numCols > numFig
    figure(numRows*numCols);
end

%New position of the figures
%Position [Left bottom width height]
screensz=get(0,'screenSize');

%screen width and height
widthSc=screensz(3);
heightSc=screensz(4);

%figure width and height
widthfig = widthSc/numCols ;
heightFig = heightSc/numRows;

figH = sort(figH);

%positions
count=0;
for j=1:numRows
    for i=1:numCols
        
        %figure
        count=count+1;
        
        %positions
        Left = (i-1)*widthfig;
        bottom = (numRows-j)*heightFig;
        scalar=0.86;
        
        %store
        Positions = [Left bottom widthfig scalar*heightFig];
        set(figH(count),'position',Positions)
    end
end



