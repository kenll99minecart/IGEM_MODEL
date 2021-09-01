function Ans = SS2(kap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
syms xm EnvZ EnvZP EnvZPR EnvZRP EnvZR OmpR OmpRP GFP RFP
kad = 0.001;
kb1 = 0.5;
kd1 = 0.5;
kpt= 1.5;
kd2 = 0.5;
kb2 = 0.5;
kph =0.05;
kd3 = 0.5;
kb3 = 0.5;
tf = 7200;

KC=20e-3;
KF=1e-3;
KF4=20e-3;
kG=0.01;
kR=0.01;
kC=0.001;
kF=0.001;
kdG=0.001;
kdR=0.001;
%equations:
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
        dxdt=[dEnvZdt;dEnvZPdt;dEnvZPRdt;dEnvZRPdt;dEnvZRdt;dOmpRdt;dOmpRPdt;dGFPdt;dRFPdt];
Ans = solve(dxdt,[EnvZ EnvZP EnvZPR EnvZRP EnvZR OmpR OmpRP GFP RFP]);
end

