
%slgetVectorStats
%
%author : steeve laquitaine
%purpose: calculate vector length and directions statistics
%         probability distribution of vector direction, each direction
%         average vector length and sem.
%
%
% nullVecDirection = [0 45 90 135 180 170 180 0];
% coor = SLpolar2cartesian(nullVecDirection',1);
% [vecdir,vecdirPdf,veclenMeanByDir,veclenSemBydir,veclenCIBydir,veclen] = slgetVectorStats(coor);

function [vecdir,vecdirPdf,veclenMeanByDir,veclenSemBydir,veclenCIBydir,veclen] = slgetVectorStats(coor)


%calculate vector lengths (norm)
veclen = sqrt(coor(:,1).^2 + coor(:,2).^2);

%calculate mean and sem of vector length by vector direction 
[~,vecdir] = SLcart2polar(coor);
vecstats = SLmakeStat(veclen,vecdir);

%re-sort vector direction and length to match
%stats outputs
vecdir = vecdir(vecstats.conditionsSubs);
veclen = veclen(vecstats.conditionsSubs);
veclenMeanByDir = vecstats.mean;
veclenSemBydir = vecstats.sem;
veclenCIBydir = cell2mat(SLmakeColumn(vecstats.confInterv));
vecdirPdf = vecstats.count/sum(vecstats.count);




