


%author: steeve laquitaine
%  date: 151014
%
%usage:
%
%        slProject('tutorial')
%
%reference: http://stackoverflow.com/questions/29187529/how-to-plot-the-projection-of-a-vector

function slProject(varargin)

%tutorial
if any(strcmp(varargin,'tutorial'))
  
  %project vector x on vector w
  x = [6 7];
  w = [1 4];
  p = (dot(x,w)/dot(x,w))*w;

  figure('color','w');
  quiver(0,0,x(1), x(2));
  hold on;
  quiver(0,0,w(1), w(2));
  quiver(0,0,p(1), p(2));
  
end