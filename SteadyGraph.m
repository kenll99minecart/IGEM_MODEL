function Varscar = SteadyGraph(kap1,kap2)
%UNTITLED2 Summary of this function goes here
%This program has been running for ... idk PLEASE ALLOCATE TIME FOR 15 MINUTES
%   Default initial
%For normal TCS
EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.17;
OmpRi=6;
initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi];

%For TCS with scaffold
EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.17;
OmpRi=6;
Sci=12;
ScEnvZPi=0;
ScEnvZPRi=0;
ScEnvZRPi=0;
ScEnvZRi=0;
ScOmpRPi=0;
ScEnvZi=0;
ScOmpRi=0;
initial2 = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;Sci;ScEnvZPi;ScEnvZPRi;ScEnvZRPi;ScEnvZRi;ScOmpRPi;ScEnvZi;ScOmpRi];

kaprange=linspace(kap1,kap2,50);
results=zeros(7,50);
resultsScar=zeros(15,50);
% %first iteration for normal
% results(:,1)=TCS(kaprange(1),0,initial);
% initial3=results(:,1);
% %for loop
% for i = 2:50
% results(:,i)=TCS(kaprange(i),0,initial3);
% initial3=results(:,i);
% 
% end
% 
% %first iteration for scaffold
% resultsScar(:,1)=TCSwscar(kaprange(1),0,initial2);
% initial4=resultsScar(:,1);
% %for loop
% for i = 2:50
% resultsScar(:,i)=TCSwscar(kaprange(i),0,initial4);
% initial4=resultsScar(:,i);
% 
% end
%20:20-20:30 10 fucking minutes
Varscar=zeros(15,10);
%for loop
initial2(8)=0;
for i = 1:10
Varscar(:,i)=TCSwscar(0.6,0,initial2);
initial2(8)=initial2(8)+1;
end
hold on;
% plot(kaprange,results(7,:),'b-');
% plot(kaprange,resultsScar(15,:),'c-');
% legend('OmpRP','OmpRPwscar');
% xlabel('kap, s/-1');
% ylabel('OmpRP/\muM');

plot(linspace(0,9,10),Varscar(7,:));
xlabel('scaf,\muM');
ylabel('OmpRP/\muM');
% resultperOMPRP=results(7,:)./OmpRi.*100;
% plot(kaprange,resultperOMPRP,'r-');
% xlabel('kap, s/-1');
% ylabel('OmpRP/\muM');
end

