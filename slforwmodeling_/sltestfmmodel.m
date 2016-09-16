

%sltestfmmodel.m
%
%author: steeve laquitaine
%
%
%inputs
%
%       fm.W_tr : is a Nv voxels x Nk channels matrix


function fm = sltestfmmodel(instances,svec_test,fm)

%calculate the channel responses using the trained weights
%and the test dataset instances
%fm.W_tr is a Nv voxels x Nk channels matrix
fm.Cresp = instances*pinv(fm.Wtrained)'; 

%re-center the channel responses relative to the match between the channel
%preference and the displayed orientation and average them
s_u = unique(svec_test);
fm.center = s_u(round(length(s_u)/2));
o = slgetmeanCresp(fm.Cresp,svec_test,fm.phi_k,fm.center);
fm.alltestCresp = o.alltestCresp;
fm.meantestCresp = o.meantestCresp;
