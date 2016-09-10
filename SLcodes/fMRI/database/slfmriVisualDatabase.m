
%slfmriVisualDatabase.m
%
%
% author: steeve laquitaine
%   date: 160729
%purpose: display matrices of voxel responses sorted by trial, clustered by 
%         voxel selectivities and another variable (subject switching to 
%         prior or evidence)
%
%         
%  usage:
%
%       slfmriInitAnalysisTaskDotDirfMRI05
%       [d,o,behbySess] = slfmriGetDBoverSessions(o,varargin);
%       [~,VoxelsParams] = slfmriGetVoxTuningParams('vonMisesprior','neuralvector',d,o);
%       slfmriVisualDatabase(d.instances(d.mySwitch==1,:),'x1',d.myRandomDir(d.mySwitch==1,:),'x2',VoxelsParams.modes1);
%       slfmriVisualDatabase(d.instances(d.mySwitch==2,:),'x1',d.myRandomDir(d.mySwitch==2,:),'x2',VoxelsParams.modes2)
%
%
%  input:
%
%       'x1', rowvar : vector variable to sort database rows
%       'x2', colvar : vector variable to sort database col

function [d_s_sel,v_s_sel] = slfmriVisualDatabase(d,varargin)

%case variable v to sort database columns
%Segmenting columns by variable v. Columns 
%with similar v are grouped together within 
%two vertical black bars
if any(strcmp(varargin,'x2')) && ~any(strcmp(varargin,'x1'))
    
    fprintf('%s \n','(slfmriVisualDatabase) Segmenting columns by variable v. Columns with similar v are grouped together within two vertical black bars')
    
    v = varargin{find(strcmp(varargin,'x2'))+1};
    
    %sort database by columns values
    [v_s,ix] = sort(v);
    d_s = d(:,ix);
    
    %isolate columns without values
    notsel = isnan(v_s);
    sel = ~isnan(v_s);
    
    %number of trials
    nTrials = size(d_s(:,sel),1);
    
    %visualize database of columns with values
    figure('color','w');
    d_s_sel = d_s(:,sel);
    imagesc(d_s_sel); box off

    %visuzalize groups of columns with similar values
    BinStart = [1:10:351];
    BinEnd = 10:10:360;   
                
    %delimitate columns with different value group
    v_s_sel = v_s(sel);
    for i = 1 : length(BinStart)        
        voxCols = SLfindBetween(v_s_sel,BinStart(i),BinEnd(i),'includeBinEnd');
        if ~isempty(voxCols)
            hold on; vline(max(voxCols),'-k')
        end
    end
    
    nCol_sel = size(d_s_sel,2);
    set(gca,'xtick',1:20:nCol_sel,'xticklabel',v_s_sel(1:20:nCol_sel))   
    xlabel('x2')
end

%case variable x1 and x2 to sort database rows and columns
if any(strcmp(varargin,'x1')) && any(strcmp(varargin,'x2'))
          
    %sort database rows by x1 values
    fprintf('%s \n','(slfmriVisualDatabase) Sorting database rows by variable x1.')    
    v1 = varargin{find(strcmp(varargin,'x1'))+1};
    [v1_s,ix] = sort(v1);
    d_s = d(ix,:);    
    
    %sort and group database by columns with similar values x2
    v2 = varargin{find(strcmp(varargin,'x2'))+1};
    [d_s_sel,v_s_sel] = slfmriVisualDatabase(d_s,'x2',v2);
        
    nRow = size(d_s_sel,1);
    set(gca,'ytick',1:1:nRow,'yticklabel',v1_s(1:1:nRow))  
    ylabel('x1')
end
