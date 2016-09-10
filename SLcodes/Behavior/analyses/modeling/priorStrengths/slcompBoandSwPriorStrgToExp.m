
%slcompBoandSwPriorStrgToExp.m
%
%
% author: steeve laquitaine
%purpose: check the null hypothesis that prior strengths for the 
%         Bayesian observer and the switching observers are same as the
%         experimental prior strengths (accurate).
%         non parametric test : p < 0.05 indicates that the null 
%         hypothesis that the strengths are the same is rejected : they are
%         the different
%         null hyp: data in the vector X come from a distribution whose 
%         median is scalar M.
%         H=1 null can be rejected.
%         
% 
%usage :
% 
%        cd('~/Dropbox/myDropbox/Codes/projInference/data/') 
%        [~,~, data] = slcsvRead('data_PriorStrngths.csv'); 
%        [Pbo,Psw,subBOstd,subswstd,subExpstd] = slIsPriorStrgDisNormal(data);
%        [pbo,psw,Hbo,Hsw,wbo,wsw] = slcompBoandSwPriorStrgToExp(subBOstd,subswstd,subExpstd);


function [pbo,psw,Hbo,Hsw,wbo,wsw] = slcompBoandSwPriorStrgToExp(subBOstd,subswstd,subExpstd)

%get unique experimental priors
Expstd = subExpstd(:,1);

%compare each experimental prior with the distribution of prior strengths
%obtained from each model.
for i = 1 : 4
    [pbo(i),Hbo(i),wbo(i)] = signrank(subBOstd(i,:),Expstd(i));
    [psw(i),Hsw(i),wsw(i)] = signrank(subswstd(i,:),Expstd(i));
end
wbo = [wbo.signedrank];
wsw = [wsw.signedrank];
