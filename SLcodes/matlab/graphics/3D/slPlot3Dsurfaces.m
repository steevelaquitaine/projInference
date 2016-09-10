

%slPlot3Dsurfaces.m
%
%
% author: steeve laquitaine
%   date: 150724
%purpose: plot surfaces. matlab surf function is not intuitive enough.
%
%
%   usage:
%
%           slPlot3Dsurfaces(1:10,1:10,rand(10,10),ones(10,10))
%           slPlot3Dsurfaces(1:8,1:10,rand(10,8),ones(10,8))


function slPlot3Dsurfaces(x,y,z,C)

%checks
if any(size(z,1)==1)
   fprintf('(slpot3Dsurfaces) I"m sorry but z should be a matrix not a vector \n'); 
end

if size(z,2) ~= length(x)
    fprintf('(slpot3Dsurfaces) z row nb must equal length of x \n'); 
end

if size(z,1) ~= length(y)
    fprintf('(slpot3Dsurfaces) z col nb must equal length of y \n');  
end

if any(size(C) ~= size(z))
     fprintf('(slpot3Dsurfaces) z and c must have same size \n');  
end

set(gcf,'color','w')
surf(x,y,z,C)
shading interp
