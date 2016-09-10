tic
load data
kli(coh==0.24)=k(1);
kli(coh==0.12)=k(2);
kli(coh==0.06)=k(3);
kpi(pstd==80)=k(4);
kpi(pstd==40)=k(5);
kpi(pstd==20)=k(6);
kpi(pstd==10)=k(7);
[logL,ML]=GirshickML(upo,d,kpi,kli);
-sum(logL)
toc