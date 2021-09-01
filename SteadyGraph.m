function resultSS = SteadyGraph(kap1,kap2,scaf1,scaf2,surf1,time)
%recommend: kap2=0.01
%This program has been running for ... idk PLEASE ALLOCATE TIME FOR 15 MINUTES
%% Clearing the memory
close all; % Closes the figures
clc; % Clears the command window
hold off

%% default initialization
%   Default initial
%For normal TCS
EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.1;
OmpRi=6;
GFPi=0;
RFPi=0;
%;GFPi;RFPi
initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi];

%For anti-TCS
EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.1;
OmpRi=6;
GFPi=0;
RFPi=0;
%;GFPi;RFPi
initials = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi];

%default initial:
initialg=zeros(1,9);

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

%default initial2:
initial2g=ones(1,17).*5;%5
%initial2g(2:5)=0;

datapoints=10;
kaprange=linspace(kap1,kap2,datapoints);
results=zeros(9,datapoints);
resultsScar=zeros(18,datapoints);
if(time)
%first iteration for normal
results(:,1)=TCS(kaprange(1),0,initial,0);
initial3=results(:,1);
%for loop
for i = 2:datapoints
results(:,i)=TCS(kaprange(i),0,initial3,0);
initial3=results(:,i);

end

%first iteration for scaffold
resultsScar(:,1)=TCSwscar(kaprange(1),0,initial2,0);
initial4=resultsScar(:,1);
%for loop
for i = 2:datapoints
resultsScar(:,i)=TCSwscar(kaprange(i),0,initial4,0);
initial4=resultsScar(:,i);

end
%20:20-20:30 10 fucking minutes
Varscar=zeros(17,10);
%for loop
initial2(8)=0;
for i = 1:10
Varscar(:,i)=TCSwscar(0.5,0,initial2,0);
initial2(8)=initial2(8)+1e-5;
end
figure(1);
plot(kaprange,results(7,:),'b-',kaprange,resultsScar(7,:),'r-');
legend('OmpRP','OmpRPwscar');
xlabel('kap, s/-1');
ylabel('OmpRP/\muM');

figure(2);
plot(linspace(0,50e-5,10),Varscar(7,:));
xlabel('scaf,\muM');
ylabel('OmpRP/\muM');
% resultperOMPRP=results(7,:)./OmpRi.*100;
% plot(kaprange,resultperOMPRP,'r-');
% xlabel('kap, s/-1');
% ylabel('OmpRP/\muM');

if(surf1)
mat=zeros(15,15);
initial5=initial2;
initial5(8)=scaf1;
resultsn=zeros(15,15);
kaprange2=linspace(kap1,kap2,15);
kapplot=linspace(kap2,kap1,15);
for i = 1:15
%for loop
for j = 1:15
resultsn(:,j)=TCSwscar(kaprange2(j),0,initial5,0);
end
mat(i,:)=resultsn(7,:);
initial5(8)=initial5(8)+((scaf2-scaf1)/14);
end
figure(4);
x=linspace(scaf1,scaf2,15);
mat

mesh(x,kapplot,mat);
ylabel('kap, s/-1');
xlabel('scaf, \muM');
zlabel('OmpRP,\muM');
end



else %if(time)
resultSS=zeros(7,101);
resultSSscaf=zeros(16,101);
resultanti=zeros(13,101);
kaprangeSS=linspace(kap1,kap2,101);

%iteration for changing intial guess for TCS
resultSS(:,1)=TCS(kaprangeSS(1),0,initial,0,initialg);
initialchange=resultSS(:,1);
%simulate once using ode45

for i=2:101
% resultSS(:,i)=TCS(kaprangeSS(i),0,initial,1,initialchange);
% initialchange=resultSS(:,i);
resultSS(:,i)=TCS(kaprangeSS(i),0,initial,0,initialchange);
initialchange=resultSS(:,i);
end
% for i=2:31
% if(mod(i,10)==0)
% resultSS(:,i)=TCS(kaprangeSS(i),0,initial,1,initialchange);
% initialchange=resultSS(:,i);
% else
% resultSS(:,i)=TCS(kaprangeSS(i),0,initial,1,initialchange);
% initialchange=resultSS(:,i);
% end
% end

% for i=1:51
% resultSS(:,i)=TCS(kaprangeSS(i),0,initialg,1,initial);
% end
resultanti(:,1)=TCSanti(kaprangeSS(1),0,initials,0,initialg);
initialchange3=resultanti(:,1);
for i=2:101
% resultSS(:,i)=TCS(kaprangeSS(i),0,initial,1,initialchange);
% initialchange=resultSS(:,i);
resultanti(:,i)=TCSanti(kaprangeSS(i),0,initials,0,initialchange3);
initialchange3=resultanti(:,i);
end
fprintf('done TCS');

%*************************************************
%iteration for changing intial guess for TCS with scarfold
resultSSscaf(:,1)=TCSwscar(kaprangeSS(1),0,initial2,0,initial2g);
%simulate once using ode45
initialchange2=resultSSscaf(:,1);
%*************************************************


% for i=2:31
% if(mod(i,10)==0)
% resultSSscaf(:,i)=TCSwscar(kaprangeSS(i),0,initial2,1,initialchange2);
% 
% else
% resultSSscaf(:,i)=TCSwscar(kaprangeSS(i),0,initial2,1,initialchange2);
% end
% 
% initialchange2=resultSSscaf(:,i);
% end

%*************************************************
for i=2:101 
resultSSscaf(:,i)=TCSwscar(kaprangeSS(i),0,initial2,0,initialchange2);
initialchange2=resultSSscaf(:,i);
count=i
end
%**************************************************

% for i=1:51
% resultSSscaf(:,i)=TCSwscar(kaprangeSS(i),0,initial2,1,initial2g);
% %TotalEnvZ=initial2(1)+initial2(2)+initial2(3)+initial2(4)+initial2(5)+initial2(9)+initial2(10)+initial2(11)+initial2(12)+initial2(13)-resultSSscaf(1,i)-resultSSscaf(2,i)-resultSSscaf(3,i)-resultSSscaf(4,i)-resultSSscaf(5,i)-resultSSscaf(9,i)-resultSSscaf(10,i)-resultSSscaf(11,i)-resultSSscaf(12,i)-resultSSscaf(13,i)      
% %fprintf('next');
% end

figure(5);
plot(kaprangeSS,resultSS(7,:),'r-',kaprangeSS,resultSSscaf(7,:),'b-',kaprangeSS,resultSSscaf(7,:)+resultSSscaf(15,:))
%resultSSscaf;
legend('OmpRP','OmpRP with scaffold','OmpRP+OmpRP.SC');
xlabel('kap, s/-1');
ylabel('OmpRP/\muM');

scafrange=linspace(scaf1,scaf2,51);
diffscaf=zeros(17,51);
teminitial=initial2;
% for i=1:51
% teminitial(8)=scafrange(i); 
% diffscaf(:,i)=TCSwscar(0.7,0,teminitial,0,initial2g);
% end

% figure(6);
% plot(scafrange,diffscaf(7,:),'b-')
% %resultSSscaf;
% legend('OmpRP with scaffold');
% xlabel('scaf, \muM');
% ylabel('OmpRP/\muM');

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

OmpRPn=resultSS(7,:);
GFPmn=kC.*OmpRPn.^2./(OmpRPn.^2+KC.^2)./kdGm;
RFPmn=kF.*OmpRPn.^2./(OmpRPn.^2+KF.^2).*(1-(OmpRPn.^2./(OmpRPn.^2+KF4.^2)))./kdRm;
antiGn=RFPmn;
antiRn=GFPmn;
leftGFPn=GFPmn.*KdG./(KdG+antiGn);
leftRFPn=RFPmn.*KdR./(KdR+antiRn);

GFPantin=kG.*leftGFPn./kdG;
RFPantin=kR.*leftRFPn./kdR;

OmpRP=resultSSscaf(7,:)+resultSSscaf(15,:);
GFPm=kC.*OmpRP.^2./(OmpRP.^2+KC.^2)./kdGm;
RFPm=kF.*OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)))./kdRm;
antiG=RFPm;
antiR=GFPm;
leftGFP=GFPm.*KdG./(KdG+antiG);
leftRFP=RFPm.*KdR./(KdR+antiR);

GFPanti=kG.*leftGFP./kdG;
RFPanti=kR.*leftRFP./kdR;

OmpRPt=resultSS(7,:);
GFPSS=kG.*kC.*OmpRPt.^2./(OmpRPt.^2+KC.^2)./kdG;
RFPSS=kR.*kF.*OmpRPt.^2./(OmpRPt.^2+KF.^2).*(1-(OmpRPt.^2./(OmpRPt.^2+KF4.^2)))./kdR;

OmpRPf=resultSSscaf(7,:)+resultSSscaf(15,:);
GFPSSF=kG.*kC.*OmpRPf.^2./(OmpRPf.^2+KC.^2)./kdG;
RFPSSF=kR.*kF.*OmpRPf.^2./(OmpRPf.^2+KF.^2).*(1-(OmpRPf.^2./(OmpRPf.^2+KF4.^2)))./kdR;

figure(7);
plot(kaprangeSS,GFPSS,kaprangeSS,RFPSS,kaprangeSS,GFPSSF,kaprangeSS,RFPSSF);
legend('GFP only TCS','RFP only TCS','GFP with scaffold','RFP with scaffold');

figure(8);
plot(kaprangeSS,GFPantin,kaprangeSS,RFPantin,kaprangeSS,GFPSS,kaprangeSS,RFPSS)
legend('GFP with antisense','RFP with antisense','GFP only TCS','RFP only TCS')

figure(9);
plot(kaprangeSS,GFPSSF,kaprangeSS,RFPSSF,kaprangeSS,GFPanti,'g-',kaprangeSS,RFPanti,'r-');
legend('GFP with scaffold','RFP with scaffold','GFP with scaffold & antisense','RFP with scaffold & antisense');

figure(10);
plot(kaprangeSS,GFPantin,kaprangeSS,RFPantin,kaprangeSS,GFPanti,kaprangeSS,RFPanti);
legend('GFP with antisense','RFP  with antisense','GFP with antisense & scaffold','RFP with antisense & scaffold');
%state switches at around 0.02\muM 
%3:06-3:25(0,0.01)
% figure(8);
% plot(kaprangeSS,resultSS(8,:),kaprangeSS,resultSS(9,:),kaprangeSS,GFPSS,kaprangeSS,RFPSS);
% legend('GFP only TCS','RFP only TCS','GFPSS only TCS','RFPSS only TCS');

figure(11);
plot(RFPSS,GFPSS,'r',RFPantin,GFPantin,'g',RFPSSF,GFPSSF,'b',RFPanti,GFPanti,'m');
legend('only TCS','only antisense','only scaffold','scaffold & antisense');
xlabel('RFP');
ylabel('GFP');

%% anti2 without scaffold
GFPanti2=zeros(101,1);
RFPanti2=zeros(101,1);

EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.1;
OmpRi=6;
GFPi=0;
RFPi=0;
initial3 = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi];

fprintf('now anti2');
for i=1:101
Tem=TCSanti2(kaprangeSS(i),0,initial3);
GFPanti2(i)=Tem(15);
RFPanti2(i)=Tem(17);
count=i
end
GFPanti2
RFPanti2
%figure(12);
%plot(kaprangeSS,GFPantin,kaprangeSS,RFPantin,kaprangeSS,GFPanti2,kaprangeSS,RFPanti2)
%legend('GFP with antisense','RFP with antisense','GFP with anti2','RFP with anti2')
hold off;
figure(15);
plot(kaprangeSS,GFPanti2,'g',kaprangeSS,RFPanti2,'r');
legend('GFP','RFP');

matGFPanti2=[kaprangeSS;GFPanti2'];
matRFPanti2=[kaprangeSS;RFPanti2'];
matGFPantin=[kaprangeSS;GFPantin];
matRFPantin=[kaprangeSS;RFPantin];
P2=InterX(matGFPanti2,matRFPanti2);
P1=InterX(matGFPantin,matRFPantin);
figure(16);
plot(kaprangeSS,GFPanti2,'g',kaprangeSS,RFPanti2,'r',kaprangeSS,GFPantin,kaprangeSS,RFPantin,P2(1,:),P2(2,:),'bo',P1(1,:),P1(2,:),'ro');
legend('GFP with anti2','RFP with anti2','GFP with anti','RFP with anti');

P3=InterX([kaprangeSS;GFPSS],[kaprangeSS;RFPSS]);
figure(17);
plot(kaprangeSS,GFPanti2,'g',kaprangeSS,RFPanti2,'r',kaprangeSS,GFPSS,kaprangeSS,RFPSS,P2(1,:),P2(2,:),'bo',P3(1,:),P3(2,:),'bo');
legend('GFP with anti2','RFP with anti2','GFP','RFP');

%% anti2 with scaffold
GFPanti2scaf=zeros(101,1);
RFPanti2scaf=zeros(101,1);

EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.1;
OmpRi=6;
GFPi=0;
RFPi=0;
initial4 = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;Sci;ScEnvZi;ScEnvZPi;ScEnvZPRi;ScEnvZRPi;ScEnvZRi;ScOmpRPi;ScOmpRi;EnvZD;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi];

fprintf('now anti2 with scaffold');
for i=1:101
Tem=TCSwscafd2(kaprangeSS(i),0,initial4);
GFPanti2scaf(i)=Tem(24);
RFPanti2scaf(i)=Tem(26);
count=i
end

hold off;
figure(18);
plot(kaprangeSS,GFPanti2scaf,'g',kaprangeSS,RFPanti2scaf,'r');
legend('GFP with anti2 & scaffold','RFP with anti2 & scaffold');

matGFPanti2scaf=[kaprangeSS;GFPanti2scaf'];
matRFPanti2scaf=[kaprangeSS;RFPanti2scaf'];
matGFPanti=[kaprangeSS;GFPanti];
matRFPanti=[kaprangeSS;RFPanti];
P4=InterX(matGFPanti2scaf,matRFPanti2scaf);
P5=InterX(matGFPanti,matRFPanti);
figure(19);
plot(kaprangeSS,GFPanti2scaf,'g',kaprangeSS,RFPanti2scaf,'r',kaprangeSS,GFPanti,kaprangeSS,RFPanti,P4(1,:),P4(2,:),'bo',P5(1,:),P5(2,:),'ro');
legend('GFP with anti2 & scaffold','RFP with anti2 & scaffold','GFP with anti & scaffold','RFP with anti & scaffold');

figure(20);
plot(kaprangeSS,GFPanti2scaf,'g',kaprangeSS,RFPanti2scaf,'r',kaprangeSS,GFPSS,kaprangeSS,RFPSS,P4(1,:),P4(2,:),'bo',P3(1,:),P3(2,:),'bo');
legend('GFP with anti2','RFP with anti2','GFP','RFP');

end %if(time)
end

