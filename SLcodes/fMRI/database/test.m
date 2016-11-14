



%%
%load data, shape and classify
i=1;
subs{i} = 's25';
roi = 'V1';
load(['~/data/datafMRI/sltaskdotdirfmri05/'...
    'slStckAnalyses/Concatenation/classif/' subs{i} '/prior225/myRandomCoh/'...
    'accAtTimeleaveOneOutfisherbalancByRemovI/' roi '/d.mat'])

%sort data into cells for each variable
[iCells,v_u,iCells_ix] = slmat2CellbyVarSortByVar(d.instances,d.mySwitch,d.myRandomCoh,0.06);
iCells = iCells(v_u==1 |  v_u==2);

%get variable 2 associated with each cell rows
iCells_ix = iCells_ix(v_u==1 |  v_u==2);
v_u = [1 2];
for i = 1 : length(v_u)
    var2{i} = d.myRandomDir(iCells_ix{i});
end

%balance variable 2 between dataset cells
[dataCells_bal,ixKeep,v_u] = slbalanceDataForVar(iCells,var2);
c = leaveOneOut(dataCells_bal);
