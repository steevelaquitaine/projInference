

% author: steeve laquitaine
%   date: 150726
%purpose: create times series with defined statistics
%
%
%  usage: 
%
%   
%       series = slMakeStatSeries([15 85 155 225 295],[18 18 41 67 41]);

function series = slMakeStatSeries(vals,rep)

series = [];

for i = 1 : length(vals)   
    
   thisSeries  = repmat(vals(i),1,rep(i));
   series = [series thisSeries];
   
end