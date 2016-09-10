%when motor noise is a von mises.
llh=vmPdfs(1:1:360,360,2.3,'norm');
Mn=vmPdfs(1:1:360,0,2.3,'norm');
cconv1=circConv(llh,Mn);
%when motor noise is flat.
llh=vmPdfs(1:1:360,360,2.3,'norm');
Mn=vmPdfs(1:1:360,0,0,'norm');
cconv2=circConv(llh,Mn);
%there is no motor noise.
llh=vmPdfs(1:1:360,360,2.3,'norm');
Mn=vmPdfs(1:1:360,0,600,'norm');
cconv3=circConv(llh,Mn);