

%slfmriEqualizeDirDist.m
%
%
% author: steeve laquitaine
%purpose: equalize direction distribution between conditions (e.g., switching)
%
%  usage: 
%
%       slfmriInitAnalysisTaskDotDirfMRI05
%       d = slsimfMRIdatabase
%       d_dirBalbySw = slfmriEqualizeDirDist(d,o)


function [d_dirBalbySw,coh_u] = slfmriEqualizeDirDist(d,o)

% 
%% Equalize distribution of directions across switching variable
du = unique(d.myRandomDir);
nDu = length(du);
coh_u = unique(d.myRandomCoh);

% for icoh = 1: length(coh_u)
%     ixDiriBal_sw1 = [];
%     ixDiriBal_sw2 = [];
%     coh_i = coh_u(icoh);
%     parfor i = 1 : nDu
%         %trials this direction & condition
% %         ixDiri_Sw1 = find(d.mySwitch==1 & d.myRandomCoh==coh_i & d.myRandomDir==du(i));
%         ixDiri_Sw2 = find(d.mySwitch==2 & d.myRandomCoh==coh_i & d.myRandomDir==du(i));
%         %equalize number of occurence for each direction
%         %by removing additional directions
%         minNb = min([length(ixDiri_Sw1); length(ixDiri_Sw2)]);
%         %we keep the first n min trials
%         ixDiriBal_sw1 = [ixDiriBal_sw1; ixDiri_Sw1(1:minNb)];
%         ixDiriBal_sw2 = [ixDiriBal_sw2; ixDiri_Sw2(1:minNb)];
%     end 
%     %re-build the database
%     d_dirBalbySw{icoh} = d;
%     fieldnm = fieldnames(d);
%     ixAll = sort([ixDiriBal_sw1;ixDiriBal_sw2]);
%     for i = 1:length(fieldnames(d))
%         d_dirBalbySw{icoh}.(fieldnm{i}) = d.(fieldnm{i})(ixAll,:);
%     end
% end


for icoh = 1: length(coh_u)
    ixDiriBal_sw1 = [];
    ixDiriBal_sw2 = [];
    coh_i = coh_u(icoh);
    parfor i = 1 : nDu
        %trials this direction & condition
%         ixDiri_Sw1 = find(d.mySwitch==1 & d.myRandomCoh==coh_i & d.myRandomDir==du(i));
        ixDiri_Sw2 = find(d.mySwitch==2 & d.myRandomCoh==coh_i & d.myRandomDir==du(i));
        %equalize number of occurence for each direction
        %by removing additional directions
        minNb = min([length(ixDiri_Sw2)]);
        %we keep the first n min trials        
        ixDiriBal_sw2 = [ixDiriBal_sw2; ixDiri_Sw2(1:minNb)];
    end 
    %re-build the database
    d_dirBalbySw{icoh} = d;
    fieldnm = fieldnames(d);
    ixAll = sort([ixDiriBal_sw2]);
    for i = 1:length(fieldnames(d))
        d_dirBalbySw{icoh}.(fieldnm{i}) = d.(fieldnm{i})(ixAll,:);
    end
end


%plot distribution of directions sorted by switching variable
% figure('color','w'); 
% for icoh = 1 : length(coh_u)    
%     subplot(1,2,icoh); hold all
%     h1 = hist(d_dirBalbySw{icoh}.myRandomDir(d_dirBalbySw{icoh}.mySwitch==1),1:1:360);
%     bar([1:1:360]-4,h1,10,'facecolor',[0 0.5 .7],'lineStyle','none')
%     h2 = hist(d_dirBalbySw{icoh}.myRandomDir(d_dirBalbySw{icoh}.mySwitch==2),1:1:360);
%     bar([1:1:360]+4,h2,10,'facecolor','k','lineStyle','none')
%     ylabel('Count'); xlabel('Motion direction (deg)'); xlim([0 360])
%     legend('Sw-p','Sw-d'); alpha(0.9)
%     set(gca,'xtick',unique(d_dirBalbySw{icoh}.myRandomDir),'xticklabel',unique(d_dirBalbySw{icoh}.myRandomDir))
%     vline(o.priormean,':k') ; title(['Direction dist. by switching for' num2str(coh_u(icoh)) '% coh']); box off   
% end



%plot distribution of directions sorted by switching variable
figure('color','w'); 
for icoh = 1 : length(coh_u)    
    h2 = hist(d_dirBalbySw{icoh}.myRandomDir(d_dirBalbySw{icoh}.mySwitch==2),1:1:360);
    bar([1:1:360]+4,h2,10,'facecolor','k','lineStyle','none')
    ylabel('Count'); xlabel('Motion direction (deg)'); xlim([0 360])
    legend('Sw-p','Sw-d'); alpha(0.9)
    set(gca,'xtick',unique(d_dirBalbySw{icoh}.myRandomDir),'xticklabel',unique(d_dirBalbySw{icoh}.myRandomDir))
    vline(o.priormean,':k') ; title(['Direction dist. by switching for' num2str(coh_u(icoh)) '% coh']); box off   
end