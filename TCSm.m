function Ans = TCSm(x)
% important:  according to literature, EnvZ must be far greater than OmpR

%x=gaminv(x,1,0.1); %assume gamma distribution
x=unifinv(x,0,1);
%parameters goes here
kap=0.1;
kad = x(1);%0.001
kb1 = x(2);
kd1 =x(3);
kpt=x(4); %1.5/0.0015 %need to sensitivity analysis
kd2 =x(5);
kb2 = x(6);
kph =x(7);%0.0015e-6%wrong prev 0.05
kd3 = x(8);
kb3 = x(9);
tf = 100000;
x
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
EnvZi=0.1;
OmpRi=6;
GFPi=0;
RFPi=0;
initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi];

opt=odeset('RelTol',1e-5,'NonNegative',1,'Events',@reachSS);%
[tt, NN] = ode15s(@f, [0; tf], initial,opt);
Ans=NN(end,7);
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
        dOmpRdt=-2.*kb1.*EnvZP.*(OmpR.^2)+2.*kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
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
      %fprintf('SS \n');
      position = (maxChange < 1e-10) - 0.5;        
      isterminal = 1;        
      direction = 1;    
    end  
end

