load data
tic
kp1=2; %checked and it's ok
kp2=0:0.1:9.9;
kp3=0:2:18;
kp4=0:2:18;
kli(coh==0.24)=k(1);
kli(coh==0.12)=k(2);
kli(coh==0.06)=k(3);
kpi(pstd==80)=kp1;
tic
for j=1:1:100
    for k=1:10
        for l=1:10
            kpi(pstd==40)=kp2(j);
            kpi(pstd==20)=kp3(k);
            kpi(pstd==10)=kp4(l);
            [logL,ML]=GirshickML(upo,d,kpi,kli);
            slogL(i,j,k)=-sum(logL);
            
            hold all
            h=plot(1,1);
            plot(h,slogL(i,j,k),'ko')
            drawnow
        end
        toc
    end
end
figure;
[x,y,z]=meshgrid(1:10,1:10,1:10);
scatter3(x(:),y(:),z(:),100,slogL(:))
