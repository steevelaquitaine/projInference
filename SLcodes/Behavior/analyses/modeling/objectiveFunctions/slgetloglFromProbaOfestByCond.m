
%slgetloglFromProbaOfestByCond.m
%
%
%  author: steeve laquitaine
% purpose: calculate logl of data from model predicted estimate 
%          distribution and conditions
%
%usage : 
%
%  db = SLMakedatabank({'sub02'},'experiment','vonMisesPrior','dataPath','~/data/dataPsychophy/proj01_priorStrength/');
%  data = db.estimatesDeg; coh = db.stimStrength; disp = db.stimFeatureDeg; pstd = db.Pstd;
%  loglpertrial = slgetloglFromProbaOfestByCond(pest,cond,data,disp,coh,pstd);


function [loglpertrial,logl] = slgetloglFromProbaOfestByCond(pest,cond,data,disp,coh,pstd)

%get probability of each data 
%given best model fit
for i = 1 : length(data)   
    cond_i = find(cond(:,1)==pstd(i) & cond(:,2)==coh(i) & cond(:,3)==disp(i));
    pest_i(i) = pest(cond_i,data(i)) ;        
end

loglpertrial = log(pest_i);
logl = sum(log(pest_i));



