


function fm = slvoxfmKFoldCVdec(instances,svec,Nf)

%divide the dataset into 5 folds
%5 folds per stimulus feature (motion direction)
%to make sure data for all stimulus feature are present
%during training. The number of trials is equalized by 
%condition to prevent biasing training toward certain
%features.
%sort instances in structure by stimulus feature
[i_struc,v_u] = slmat2structByVar(instances,svec);

%balance instances between simulus features
instancesBal = setNumInstancesByClass(i_struc,'balanceByRemovI');

%get vectors of associated stimuli features
v_u = SLmakeColumn(v_u);
v_u_all = v_u(:,ones(1,size(instancesBal{1},1)))';

%get row starting and ending indices for each fold 
%by stimulus feature (col)
[f_st,f_end,~,Nivi] = slgetkFoldrows(instancesBal,Nf);

%decode from each train/test fold combinationW_tr_f(:,:,f) = W_tr;
Nv = size(instances,2);
Ni = size(instances,1);
Nk = 5;
W_tr_f = nan(Nv,Nk,Nf);
svecTest_f = nan(Nivi,Nf);
allTestCresp_t_f = nan(Nivi,Nk,Nf);

%Train and test each fold combination
for f = 1 : Nf
    
    %this fold test and train instances
    %picked up from 1/Nf of each stimulus feature
    testInst = []; trainInst = []; svecTrain = []; svecTest = [];
    for s_u = 1 : length(unique(svec))        
        
        %this fold test instances
        testInst = [testInst; instancesBal{s_u}(f_st(f) : f_end(f),:)];        
        %this fold's train instance positions
        trainIx = setdiff(1:Nivi,f_st(f) : f_end(f));
        %this fold's train instances
        trainInst = [trainInst; instancesBal{s_u}(trainIx ,:)];        
        %this fold's train stimulus feature
        svecTrain = [svecTrain; v_u_all(trainIx,s_u)];
        %this fold's test stimulus feature
        svecTest = [svecTest; v_u_all(f_st(f):f_end(f),s_u)];
    end
   
    %decode this test fold (train then test)
    fm_f = sltrainfmmodel(trainInst,svecTrain);   
    fm_f = sltestfmmodel(testInst,svecTest,fm_f);
        
    %backup this fold
    fm.Wtrained_f(:,:,f) = fm_f.Wtrained;    
    fm.svecTest_f(:,f) = svecTest;
    fm.meanTestCresp_f(f,:) = fm_f.meantestCresp;
    fm.allTestCresp_t_f(:,:,f) = fm_f.alltestCresp;
    fm.center = fm_f.center;
    
    %backup model    
    fm.s = fm_f.s;    
    fm.phi_k = fm_f.phi_k;
    fm.K = fm_f.K;
    fm.phi_k = fm_f.phi_k;
    fm.s = fm_f.s;   
    fm.f_k_s = fm_f.f_k_s;
    
end

%reshape as a matrix of Nf folds x Nivi instances rows x Nk channels
%and average channel responses over rows
fm.allTestCresp_t = reshape(permute(fm.allTestCresp_t_f,[1 3 2]),Nf*Nivi,Nk); 

%LLH mean and std
fm.meanLLH_all = mean(fm.allTestCresp_t,1);
fm.stdLLH_all = std(fm.allTestCresp_t);

%likelihood with sem over direction repeats
figure('color','w')
errorarea(fm.allTestCresp_t','k',[.5 .5 .5])
plot(fm.meanLLH_all,'k.','markersize',40)
set(gca,'xtick',1:length(v_u),'xticklabel',v_u - fm_f.center)
ylabel('Channels responses (au)')
xlabel('Channels by direction preference distance to the displayed direction (deg)')
box off
title([num2str(Nf) ' folds cross-validated channel responses (averaged over) ' num2str(Nf) 'x' num2str(Nivi) ' instances'])



