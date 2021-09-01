function [] = SteadyScaf(kap1,kap2,scafd1,scafd2)
%SS for varying degradation rate

%For TCS with scaffold
EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.1;
OmpRi=6;
Sci=0;%12
ScEnvZPi=0;
ScEnvZPRi=0;
ScEnvZRPi=0;
ScEnvZRi=0;
ScOmpRPi=0;
ScEnvZi=0;
ScOmpRi=0;
GFPi=0;
RFPi=0;
EnvZD=0;
%;GFPi;RFPi
initial2 = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;Sci;ScEnvZi;ScEnvZPi;ScEnvZPRi;ScEnvZRPi;ScEnvZRi;ScOmpRPi;ScOmpRi;EnvZD];


KC=20e-3;
KF=1e-3;
KF4=20e-3;
kG=1;
kR=1;
kC=0.01;
kF=0.01;
kdG=0.001;
kdR=0.001;
kdGm=1;
kdRm=1;
KdG=1e-2;
KdR=1e-2;

resultSSscaf=zeros(16,101);
kaprangeSS=linspace(kap1,kap2,101);
Scafdrange=linspace(scafd1,scafd2,51);
GFPanti=zeros(101,51);
RFPanti=zeros(101,51);
GFPSSF=zeros(101,51);
RFPSSF=zeros(101,51);

for j=1:51
%*************************************************
%iteration for changing intial guess for TCS with scarfold
resultSSscaf(:,1)=TCSwscafd(kaprangeSS(1),0,initial2,Scafdrange(j));
%simulate once using ode45
initialchange2=resultSSscaf(:,1);
%*************************************************
%*************************************************
for i=2:101 
resultSSscaf(:,i)=TCSwscafd(kaprangeSS(i),0,initial2,Scafdrange(j));
initialchange2=resultSSscaf(:,i);
count=i+(j-1)*101
end
%**************************************************

OmpRP=resultSSscaf(7,:)+resultSSscaf(15,:);
GFPm=kC.*OmpRP.^2./(OmpRP.^2+KC.^2)./kdGm;
RFPm=kF.*OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)))./kdRm;
antiG=RFPm;
antiR=GFPm;
leftGFP=GFPm.*KdG./(KdG+antiG);
leftRFP=RFPm.*KdR./(KdR+antiR);

GFPanti(:,j)=kG.*leftGFP./kdG;
RFPanti(:,j)=kR.*leftRFP./kdR;

OmpRPf=resultSSscaf(7,:)+resultSSscaf(15,:);
GFPSSF(:,j)=kG.*kC.*OmpRPf.^2./(OmpRPf.^2+KC.^2)./kdG;
RFPSSF(:,j)=kR.*kF.*OmpRPf.^2./(OmpRPf.^2+KF.^2).*(1-(OmpRPf.^2./(OmpRPf.^2+KF4.^2)))./kdR;

end
figure(1);
plot(kaprangeSS,GFPSSF(:,1),'g-',kaprangeSS,RFPSSF(:,1),'r-');
figure(2);
plot(kaprangeSS,GFPSSF(:,21),'g-',kaprangeSS,RFPSSF(:,21),'r-');
figure(3);
%F(21) = struct('cdata',[],'colormap',[]);
vid=VideoWriter('v4.mp4','MPEG-4');
vid.FrameRate=10;
open(vid)
for i=1:41
axis equal
plot(kaprangeSS,GFPSSF(:,i),'g-',kaprangeSS,RFPSSF(:,i),'r-',kaprangeSS,GFPanti(:,i),kaprangeSS,RFPanti(:,i));
legend('GFP','RFP','GFP with antisense','RFP with antisense')
drawnow limitrate
F(i) = getframe(gcf);
writeVideo(vid,F(i));
end
close(vid);
movie(F);
end

