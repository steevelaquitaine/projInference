
%slfMRIplotR2density.m
%

%usage:
%
%		slfMRIplotR2density('lV1')


function slfMRIplotR2density(myRoi)


% check arguments
if ~any(nargin == [0])
 help slfMRIplotR2density
 return
end

vals1 = mlrGetOverlayValues(myRoi,'Averages',3,'pRF_left_occipital_recenter_B3','r2');
vals2 = mlrGetOverlayValues('lV1','Averages',6,'pRF_left_occipital_recenter_BT3','r2');

plot(vals1,vals2,'ko','MarkerFaceColor','k');
maxp = max([vals1; vals2]); 
hold on; plot([0 maxp],[0 maxp],'r:')
set(gcf,'color','w','position',[655 213 769 580])  
box off

title('s0311 r2: B3 vs BT3')
xlabel('Averaged Bars 3 r2');
ylabel('Averaged Bars Task 3 r2');

dison [12:59 PM]
The script is called mlrComparePRF.m

