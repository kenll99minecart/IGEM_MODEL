function Ans = TCS(kap,showplot,initial,SS,guess)
% important:  according to literature, EnvZ must be far greater than OmpR

%parameters goes here
kad = 0.0001;%0.001
kb1 = 0.5;
kd1 = 0.5;
kpt= 0.0015; %1.5
kd2 = 0.5;
kb2 = 0.5;
kph =0.05;%0.0015e-6
kd3 = 0.5;
kb3 = 0.5;
tf = 8000;

KC=20e-3;
KF=1e-3;
KF4=20e-3;
kG=1;
kR=1;
kC=0.01;
kF=0.01;
kdG=0.001;
kdR=0.001;
%kap = 0.0033.*Na/(Na+0.1463)
%EnvZPR= EnvZP.OmpR etc.
EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.17;
OmpRi=6;
GFPi=0.1;
RFPi=0.1;
%initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;GFPi;RFPi];
if(~SS)
opt=odeset('RelTol',1e-6,'NonNegative',1,'Events',@reachSS);%
[tt, NN] = ode45(@f, [0; tf], initial,opt);

Ans=NN(end,:);
else
%guess=ones(1,9).*7;
ops=optimoptions(@fsolve,'FunctionTolerance',1e-20);
ops2=optimoptions(@lsqnonlin,'FunctionTolerance',1e-20,'MaxIterations',1e10,'MaxFunctionEvaluations',1e10,'StepTolerance',1e-30,'Algorithm','trust-region-reflective');
t=1000;
%Ans=fsolve(@(x)fSS(t,x),initial,ops);
Ans=lsqnonlin(@(x)fSS(t,x),guess,zeros(1,9),[],ops2);
end

if(showplot)
hold on;
figure(1);
plot(tt,NN(:,1),tt,NN(:,2),tt,NN(:,3),tt,NN(:,4),tt,NN(:,5),tt,NN(:,6),tt,NN(:,7));
xlabel('Time / s');
ylabel('Number of moles/\muM');
legend('EnvZ','EnvZP','EnvZP.OmpRP','EnvZ.OmpRP','EnvZ.OmpR','OmpR','OmpRP');
grid on;
hold off;
% figure(2);
% plot(tt,NN(:,8),'g-',tt,NN(:,9),'r-')
% xlabel('Time / s');
% ylabel('Number of moles/\muM');
% legend('GFP','RFP');
end
    function dxdt=f(t,x)
        EnvZ=x(1);
        EnvZP=x(2);
        EnvZPR=x(3);
        EnvZRP=x(4);
        EnvZR=x(5);
        OmpR=x(6);
        OmpRP=x(7);
        %GFP=x(8);
        %RFP=x(9);
        dEnvZdt=-kap.*EnvZ+kad.*EnvZP+kd2.*EnvZRP-kb2.*EnvZ.*OmpRP+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dEnvZPdt=kap.*EnvZ-kad.*EnvZP-kb1.*EnvZP.*(OmpR.^2)+kd1.*EnvZPR;
        dEnvZPRdt=kb1.*EnvZP.*(OmpR.^2)-kd1.*EnvZPR-kpt.*EnvZPR;
        dEnvZRPdt=kpt.*EnvZPR-kd2.*EnvZRP+kb2.*EnvZ.*OmpRP-kph.*EnvZRP;
        dEnvZRdt=kph.*EnvZRP+kb3.*EnvZ.*OmpR-kd3.*EnvZR;
        dOmpRdt=-kb1.*EnvZP.*(OmpR.^2)+kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dOmpRPdt=kd2.*EnvZRP-kb2.*EnvZ.*OmpRP;
        %assume constant plasmid number CN
        %dGFPdt=kG.*kC.*OmpRP.^2./(OmpRP.^2+KC.^2)-kdG.*GFP;
        %dRFPdt=kR.*kF.*OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)))-kdR.*RFP;
        %;dGFPdt;dRFPdt
        dxdt=[dEnvZdt;dEnvZPdt;dEnvZPRdt;dEnvZRPdt;dEnvZRdt;dOmpRdt;dOmpRPdt];
    end

    function [position, isterminal, direction] = reachSS(t, X)
      Changee=abs(f(t, X));
      maxChange = max(Changee(1:7));    
      position = (maxChange < 1e-7) - 0.5;        
      isterminal = 1;        
      direction = 1;    
    end  

function dxdt=fSS(t,x)
        EnvZ=x(1);
        EnvZP=x(2);
        EnvZPR=x(3);
        EnvZRP=x(4);
        EnvZR=x(5);
        OmpR=x(6);
        OmpRP=x(7);
        GFP=x(8);
        RFP=x(9);
        dEnvZdt=-kap.*EnvZ+kad.*EnvZP+kd2.*EnvZRP-kb2.*EnvZ.*OmpRP+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dEnvZPdt=kap.*EnvZ-kad.*EnvZP-kb1.*EnvZP.*OmpR+kd1.*EnvZPR;
        dEnvZPRdt=kb1.*EnvZP.*OmpR-kd1.*EnvZPR-kpt.*EnvZPR;
        dEnvZRPdt=kpt.*EnvZPR-kd2.*EnvZRP+kb2.*EnvZ.*OmpRP-kph.*EnvZRP;
        dEnvZRdt=kph.*EnvZRP+kb3.*EnvZ.*OmpR-kd3.*EnvZR;
        dOmpRdt=-kb1.*EnvZP.*OmpR+kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dOmpRPdt=kd2.*EnvZRP-kb2.*EnvZ.*OmpRP;
        %assume constant plasmid number CN
        dGFPdt=kG.*kC.*OmpRP.^2./(OmpRP.^2+KC.^2)-kdG.*GFP;
        dRFPdt=kR.*kF.*OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)))-kdR.*RFP;
        TotalEnvZ=initial(1)+initial(2)+initial(3)+initial(4)+initial(5)-EnvZ-EnvZP-EnvZPR-EnvZRP-EnvZR;
        TotalOmpR=initial(3)+initial(4)+initial(5)+initial(6)+initial(7)-OmpR-OmpRP-EnvZPR-EnvZRP-EnvZR;
        dxdt=[dEnvZdt;dEnvZPdt;dEnvZPRdt;dEnvZRPdt;dEnvZRdt;dOmpRdt;dOmpRPdt;dGFPdt;dRFPdt;TotalEnvZ;TotalOmpR;dGFPdt;dRFPdt];
    end
end

