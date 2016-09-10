



% OUTPUT
    % equi_point
    % equi_point_opt
    % equi_point_rewardbased
    
    % plots of raw data (weighted choices)
    % plots of binned data
    % plots of fitted data (sigmoids)
    



% Find "nbHP>nbLP" & "nbHP<nbLP" in test trials (C=2)
diffBias{C,session}=[data_Bkp_aftl{Conditioni,diffB_col}]'; 
diffB_sign{C,session}=diffBias{C,session}>0;
sign=[1 0]; % positive/negative

% Prepare relevant variables
HPtmp=[];
[tmp,HPtmp]=max([RP_lft RP_rgt]');
HPlocation{C,session}=HPtmp;

% Observe choices of HP (for supp analyses)
y{C,session+2}=action==HPlocation{C,session}'+0;
    
for signidx=1:2
    
    % prepare data
    signTrials=[];
    R_set=[];
    RH=[];
    RL=[];
    VH=[];
    VL=[];
    xy=[];
    numx_unit_for_plot=[];
    PHPs_bootsp=[];
    y_binned_boottmp=[];
    R_lft=[];
    R_rgt=[];
    V_lft=[];
    V_rgt=[];
    HVlocation=[];
    xy_opt=[];
    HRlocation=[];
    xy_rewardbased=[];
    
    % Sort "nbH>nbL" & "nbH<nbL"
    signTrials=diffB_sign{C,session}==sign(signidx);
    
    % Observe choices of HP=f(number diff)
    choiceHP{signidx}=action(signTrials)==HPlocation{C,session}(signTrials)'+0;

    % Compute diffV=V(H)-V(L)
    R_set=cell2mat(data_Bkp_aftl(Conditioni,[R_lft_col R_rgt_col]));
    R_set=R_set(signTrials,:);
    clear i
    for i=1:size(R_set,1);  RH{C,session}(i)=R_set(i,HPlocation{C,session}(i));   end
    for i=1:size(R_set,1);  RL{C,session}(i)=R_set(i,3-HPlocation{C,session}(i)); end
    clear i
    VH{C,session}=HP*RH{C,session}';
    VL{C,session}=LP*RL{C,session}';
    diffV{signidx}=VH{C,session}-VL{C,session};
    
    % Compute diffB=nbH-nbL
    diffnbH_nbL{signidx}=diffBias{C,session}(signTrials);
    
    % Create plot to visualize raw data
    xy=sortrows([diffV{signidx} choiceHP{signidx} diffnbH_nbL{signidx}],1);
    x{C,session}{signidx}=xy(:,1);
    y{C,session}{signidx}=xy(:,2);
    x_diffnbH_nbL{C,session}{signidx}=xy(:,3);

    % Weight raw data (for visualization)
    xu{C,session}{signidx}=unique(x{C,session}{signidx});
    numxu=numel(xu{C,session}{signidx});
    eta=0.01;
    clear i
    for i=1:numxu
        Markr_size_y1{C,session}{signidx}(i)=sum(y{C,session}{signidx}(x{C,session}{signidx}==xu{C,session}{signidx}(i))==1)*3+eta;
        Markr_size_y0{C,session}{signidx}(i)=sum(y{C,session}{signidx}(x{C,session}{signidx}==xu{C,session}{signidx}(i))==0)*3+eta;
    end

    % bin raw data (for visualization)
    % create interval
    binwd=20;
    intv_0=-60:binwd:80;
    intv_end=intv_0(2:end);
    intv=[intv_0(1:end-1)' intv_end'];
    x_unit_binned_for_plot{C,session}=nanmean(intv');
    numx_unit_for_plot=numel(x_unit_binned_for_plot{C,session});
    clear i
    for i=1:numx_unit_for_plot;
        % find choices of HP in intervals
        HP_in_intv{C,session}{i,signidx}=y{C,session}{signidx}(x{C,session}{signidx}>=intv(i,1) & x{C,session}{signidx}<intv(i,2));
        % average choices of HP
        PHP{C,session}{i,signidx}=nanmean(y{C,session}{signidx}(x{C,session}{signidx}>=intv(i,1) & x{C,session}{signidx}<intv(i,2)));
    end
    clear i


    %%% Bootstrap data within bins to obtain std
    numboot_for_binned=1000;
    y_binned{C,session}(:,signidx)=HP_in_intv{C,session}(:,signidx);
    for j=1:numboot_for_binned % create samples
        % Consider each interval of x
        for k=1:numx_unit_for_plot
            % in case choices are observed
            if isempty(y_binned{C,session}{k,signidx})==0
                % bootstrap over choices within interval
                for i=1:numel(y_binned{C,session}{k,signidx})
                    jit=ceil(rand*numel(y_binned{C,session}{k,signidx}));
                    y_binned_boottmp{j,k}(i)=y_binned{C,session}{k,signidx}(jit);
                end
                PHPs_bootsp(j,k)=mean(y_binned_boottmp{j,k});
            else % no choices are observed
                y_binned_boottmp{j,k}=nan;
                PHPs_bootsp(j,k)=nan;
            end
        end
    end
    % compute mean and sem over the bootstrapped samples
    PHP_bootsp_for_plot{C,session}(:,signidx)=mean(PHPs_bootsp)';
    ePHP_bootsp_for_plot{C,session}(:,signidx)=std(PHPs_bootsp)';


    % Model raw choices with a sigmoid function
    x_fit=-100:0.1:100; % continuous sigmoids
    [b{C,session}(:,signidx),dev{C,session}(:,signidx),stats{C,session}(:,signidx)] = glmfit(x{C,session}{signidx},y{C,session}{signidx},'binomial','link','logit');
    yfit{C,session}(:,signidx) = glmval(b{C,session}(:,signidx), x_fit,'logit');
    equi_point{C,session}(:,signidx)=-b{C,session}(1,signidx)/b{C,session}(2,signidx);

    %%% Simulation of choices of HP for optimality
    R_lft=[data_Bkp_aftl{Conditioni,R_lft_col}]';
    R_rgt=[data_Bkp_aftl{Conditioni,R_rgt_col}]';
    V_lft=RP_lft.*R_lft;
    V_rgt=RP_rgt.*R_rgt;
    [tmp,HPlocation{C,session}]=max([RP_lft RP_rgt]');
    [tmp,HVlocation]=max([V_lft V_rgt]');
    choiceHP_opt{signidx}=[HPlocation{C,session}(signTrials)==HVlocation(signTrials)]';
    % Compute optimal equivalence point
    xy_opt=sortrows([diffV{signidx} choiceHP_opt{signidx}],1);
    x_opt{C,session}{signidx}=xy_opt(:,1);
    y_opt{C,session}{signidx}=xy_opt(:,2);
    b_opt{C,session}{signidx} = glmfit(x_opt{C,session}{signidx},y_opt{C,session}{signidx},'binomial','link','logit');
    equi_point_opt{C,session}{signidx}=-b_opt{C,session}{signidx}(1)/b_opt{C,session}{signidx}(2);


    %%% Simulation of choices of HP for choices based on reward magnitudes only.
    [tmp,HRlocation]=max([R_lft R_rgt]');
    choiceHP_rewardbased{signidx}=[HPlocation{C,session}(signTrials)==HRlocation(signTrials)]';
    % Compute equivalence point of reward-based choices (can be calculated numerically because
    % rewards are coupled)
    xy_rewardbased=sortrows([diffV{signidx} choiceHP_rewardbased{signidx}],1);
    x_rewardbased{C,session}{signidx}=xy_rewardbased(:,1);
    y_rewardbased{C,session}{signidx}=xy_rewardbased(:,2);
    b_rewardbased{C,session}{signidx} = glmfit(x_rewardbased{C,session}{signidx},y_rewardbased{C,session}{signidx},'binomial','link','logit');
    equi_point_rewardbased{C,session}{signidx}=-b_rewardbased{C,session}{signidx}(1)/b_rewardbased{C,session}{signidx}(2);
end