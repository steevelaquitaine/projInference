
% Steeve laquitaine 19032009

% need log files (with their actual name)in the dir
% rename each file as : f1.......fn
% give for each log files: 
    % p_a : smoothed probability of A1 for window of bin trials
    % A: local action


clear all
r = dir;
r(1:2,:) =[];

for i = 1:size(r); % select the files
        str = r(i).name;
        Name{i}=r(i).name;
        
       
        date{i}= str(1:8);
        proba{i}=str(10:13);
        
        
        fid = fopen(r(i).name,'r'); % id the files
        f = fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ',[18 inf]);
        f=f';

        eval(['f' int2str(i) '=f';]);
        
        %success trials
        t=find(f(:,2)==1); %trial number
        
        %data
        data=f(t,:); % success trials
        eval(['data' int2str(i) '=data';]);
        clear t

        
        A=data(:,4);   % action
        R=data(:,3);   % reward
        num_eta=1000;   
        etas=(0:num_eta-1)/num_eta;
        log_likelihood=zeros(1,num_eta);
        
        
        
        % computation of observed learning  
        sigma=5;
        range_filter=length(A);
        filter_tmp=(-range_filter:range_filter)';
        filter=exp(-filter_tmp.^2/(2*sigma^2));
        p_est=-ones(1,length(A));
        for i1=1:length(A);
            p_est(i1)=sum(A.*filter(range_filter-i1+2:2*range_filter-i1+1))./sum(filter(range_filter-i1+2:2*range_filter-i1+1));
        end
        p_est=p_est';
       
        
        
        % computation of best etas and max loglikelihood
        for i1=1:num_eta;
            p=-ones(length(A),1);
            p(1)=0.5;
            eta=etas(i1);
            for i2=1:length(A)-1;
                 p(i2+1)=p(i2)+eta*R(i2)*(A(i2)-p(i2)); 
            end
            log_likelihood(i1)=mean(A.*log(p)+(1-A).*log(1-p));    
        end
        i1=find(log_likelihood==max(log_likelihood));
        eta_max=etas(i1)
        max_log_likelihood=log_likelihood(i1)

        table=[eta_max;max_log_likelihood]






        % computation of LRIA learning 
        for i2=1:length(A)-1;
           p(i2+1)=p(i2) + eta_max*R(i2)*(A(i2)-p(i2));
        end 
        p(1)=0.5;  
        p(2:length(A),1)=p(1:end-1,1);
        
        
        
        ff=p; 
        ff(length(ff)+1:900)=NaN; % square matrix
        eval(['ff' int2str(i) ' =ff;' ]);

        fff=p_est; 
        fff(length(fff)+1:900)=NaN; % square matrix
        eval(['fff' int2str(i) ' =fff;' ]);
end  

        tablep_est = [fff1 fff2 fff3 fff4 fff5 fff6 fff7 fff8 fff9 fff10 fff11 fff12 fff13 fff14 fff15]%fff16];% fff17 fff18 fff19 fff20 fff21 fff22 fff23 fff24];% fff25 fff26 fff27 fff28 fff29 fff30 fff31 fff32 fff33 fff34 fff35 fff36 fff37 fff38] %fff39 fff40 fff41 fff42 fff43 fff44 fff45 fff46] %fff47];
        tablep = [ff1 ff2 ff3 ff4 ff5 ff6 ff7 ff8 ff9 ff10 ff11 ff12 ff13 ff14 ff15]%ff16];%ff17 ff18 ff19 ff20 ff21 ff22 ff23 ff24];% ff25 ff26 ff27 ff28 ff29 ff30 ff31 ff32 ff33 ff34 ff35 ff36 ff37 ff38] %ff39 ff40 ff41 ff42 ff43 ff44 ff45 ff46] %ff47];










