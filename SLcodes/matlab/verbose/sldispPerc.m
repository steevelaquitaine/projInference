
%sldispPerc.m
%
%
% author: steeve laquitaine
%   date: 151113
%purpose: display progress in a loop in percent
%
%  usage: put in a loop
%
%         myfun: name of your function
%          maxi: max iteration
%             i: current iteration  
% 
%      maxi = 20;
%      for i = 1 : maxi; 
%          sldispPerc('slFun',i,maxi); 
%      end

function sldispPerc(myfun,i,maxi)

%Init
if i==1
    fprintf('%s',['(' myfun ') Progress: ']);
end

%percent of the loop
p = round(i/maxi*100);
fprintf('\b\b\b\b\t%d%s',p,'%');

%end
if i == maxi
    fprintf('\n')
end
