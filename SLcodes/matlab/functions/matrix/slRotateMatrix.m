

%slRotateMatrix.m


%usage
%
%		rotMat = slRotateMatrix(myMat,90)

%in progress....


function rotMat = slRotateMatrix(myMat,th)

keyboard
myccwRotOnX = [1 　　　0 　　　　0; 
　　　　　      0 cos(th) -sin(th);
　　　　　      0 sin(th) cos(th)];


rotMat = myccwRotOnX*myMat;


