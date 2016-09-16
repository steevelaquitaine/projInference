
%slfmPlotWeights.m
%
%
% author: steeve laquitaine
%purpose: plot trained weights against true weights
%
%


function slfmPlotWeights(W_tr,W,phi_k)

figure('color','w'); 

%scatter weights
subplot(131); 
scatter(W(:), W_tr(:))
xlabel('True weights')
ylabel('Trained weights')

%trained weights
subplot(132); 
imagesc(W_tr);title('Trained weights'); 
xlabel('channels by direction preference (deg)')
ylabel('Voxels')
set(gca,'xtick',1:length(phi_k),'xticklabel',phi_k)

%true weights
subplot(133);imagesc(W); title('True weights');
xlabel('channels by direction preference (deg)')
ylabel('Voxels')
set(gca,'xtick',1:length(phi_k),'xticklabel',phi_k)


%scale to same color legend
caxis([min([W_tr(:);W(:)]) max([W_tr(:);W(:)])])