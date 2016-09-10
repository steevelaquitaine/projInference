

%slgetSameSuccessiveRuns.m
%
%author: steeve laquitaine
%usage: 
%
%       [r_cont,pr_All,pr_id,nr_id] = slgetSameSuccessiveRuns(priors,runs);
%
%outputs:
%
%       prbyPr:  the first run of all successive runs by prior (cells)
%                others runs are NaN
%       nrbyPr : same for the second run (= pr_id + 1)
%       pr_All : list of first run of all successive runs with same prior
%                non sorted by prior
%       r_cont : index vector that id contiguous runs with 0

function [r_cont,pr_All,prbyPr,nrbyPr] = slgetSameSuccessiveRuns(priors,runs)

%prior unique values
p_u = unique(priors);

%loop over priors
pr_All = [];
for i = 1 : length(p_u)    
    %get all runs for this prior
    runid_i = runs(priors==p_u(i));        
    %isolate contiguous runs
    ct = 0;
    for j = 1 : length(runid_i)-1
        if runid_i(j+1)==runid_i(j)+1
            prbyPr{i}(j) = runid_i(j);
            nrbyPr{i}(j) = runid_i(j)+1;
        else
            prbyPr{i}(j) = NaN;
            nrbyPr{i}(j) = NaN;
        end
    end
    %backup
    pr_All = [pr_All; prbyPr{i}'];    
end
%remove nan
pr_All = pr_All(~isnan(pr_All));

%create index vector to id contiguous runs with 0
r_cont = ones(size(runs));
for i = 1 : length(pr_All)
    r_cont(runs == pr_All(i)) = 0;       
    r_cont(runs == pr_All(i)+1) = 0;       
end

    
    
    
    
    
    
    