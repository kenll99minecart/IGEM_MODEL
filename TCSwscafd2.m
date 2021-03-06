function Ans = TCSwscafd2(kap,showplot,initial)
% important:  according to literature, EnvZ must be far greater than OmpR

%parameters goes here
kad = 0.0001;%0.001
kb1 = 0.5;
kd1 = 0.5;
kpt= 0.0015; %1.5/0.0015 %need to sensitivity analysis
kd2 = 0.5;
kb2 = 0.5;
kph =0.025;%0.0015e-6
kd3 = 0.5;
kb3 = 0.5;

kdScaf=0.1;
tf = 1000000;

kap2=kap;
kpt2=102.1; % ignore backwards cognate HK-RR phosphorylation (no phosphate transfer backwards)
kbLZ=0.1;
kdLZ=0.1;
kbSH3=0.4;
kdSH3=0.01;

%kap = 0.0033.*Na/(Na+0.1463)
%EnvZPR= EnvZP.OmpR etc.
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
RFPi=0;
GFPi=0;
% initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;Sci;ScEnvZPi;ScEnvZPRi;ScEnvZRPi;ScEnvZRi;ScOmpRPi;ScEnvZi;ScOmpRi;GFPi;RFPi];
%'RelTol',1e-10,

KC=20e-3;
KF=1e-3;
KF4=20e-3;
kG=0.01;%*60
kR=1;%*60
kC=0.01;
kF=0.01;
kdG=0.001;%*60
kdR=0.001;%*60


kDG=1e-3;
kDR=1e-3;
kmGtx=1e-2;
kmRtx=1e-2;
kasGtx=1e-2;
kasRtx=1e-2;
kasGb=1;
kasGd=1e-2;
kDasG=1;
kDmG=1e-2;
kasRb=1;
kasRd=1e-2;
kDasR=1;
kDmR=1e-2;
kfG=1e-2;
kfR=1e-2;
kGtl=1e-2;
kRtl=1e-2;
fdE=1e-3;
fdO=1e-3;
fpE=1e-4;
fpO=1e-4;
%kdScaf=1;
% if(~SS)
opt=odeset('RelTol',1e-5,'Events',@reachSS,'NonNegative',1);
[tt, NN] = ode15s(@f,[0; tf], initial,opt);
Ans=NN(end,:);
% else
% %guess=ones(1,17).*7;
% ops=optimoptions(@fsolve,'FunctionTolerance',1e-10);
% ops2=optimoptions(@lsqnonlin,'FunctionTolerance',1e-20,'MaxIterations',1e10,'MaxFunctionEvaluations',1e10,'StepTolerance',1e-30,'Algorithm','trust-region-reflective');
% t=1000;
% %Ans=fsolve(@(x) fSSscaf(t,x),initial,ops);,initial
% [Ans,res]=lsqnonlin(@(x)fSSscaf(t,x),guess,zeros(1,17),[],ops2);
% end


if(showplot)
hold on;
figure(1);
plot(tt,NN(:,1),tt,NN(:,2),tt,NN(:,3),tt,NN(:,4),tt,NN(:,5),tt,NN(:,6),tt,NN(:,7),tt,NN(:,8),tt,NN(:,9),tt,NN(:,10),tt,NN(:,11),tt,NN(:,12),tt,NN(:,13),tt,NN(:,14),tt,NN(:,15));
xlabel('Time / s');
ylabel('Number of moles/\muM');
 legend('EnvZ','EnvZP','EnvZP.OmpRP','EnvZ.OmpRP','EnvZ.OmpR','OmpR','OmpRP','Sc','ScEnvZ','ScEnvZP','ScEnvZPR','ScEnvZRP','ScEnvZR','ScOmpR','ScOmpRP');
grid on;
figure(2);
plot(tt,NN(:,15));
xlabel('Time / s');
ylabel('Number of moles/\muM');
legend('ScOmpRP');
figure(3);
plot(tt,NN(:,7));
xlabel('Time / s');
ylabel('Number of moles/\muM');
legend('OmpRP');
figure(4);
plot(tt,NN(:,6));
xlabel('Time / s');
ylabel('Number of moles/\muM');
legend('OmpR');
figure(5);
plot(tt,NN(:,7)+NN(:,15));
xlabel('Time / s');
ylabel('Number of moles/\muM');
legend('OmpRP+ Sc.OmpRP');


end
    function dxdt=f(t,x)
        EnvZ=x(1);
        EnvZP=x(2);
        EnvZPR=x(3);
        EnvZRP=x(4);
        EnvZR=x(5);
        OmpR=x(6);
        OmpRP=x(7);
        Sc=x(8);
        ScEnvZ=x(9);
        ScEnvZP=x(10);
        ScEnvZPR=x(11);
        ScEnvZRP=x(12);
        ScEnvZR=x(13);
        ScOmpR=x(14);
        ScOmpRP=x(15);
        EnvZD=x(16);
        mGFP=x(17);
        asGFP=x(18);
        asmGFP=x(19);
        mRFP=x(20);
        asRFP=x(21);
        asmRFP=x(22);
        GFPp=x(23);
        GFPa=x(24);
        RFPp=x(25);
        RFPa=x(26);
        %GFP=x(16);
        %RFP=x(17);
        %TCS system
        dEnvZdt=-kap.*EnvZ+kad.*EnvZP+kd2.*EnvZRP-kb2.*EnvZ.*OmpRP+kd3.*EnvZR-kb3.*EnvZ.*OmpR-kbSH3.*Sc.*EnvZ+kdSH3.*ScEnvZ-kbSH3.*ScOmpR.*EnvZ+kdSH3.*ScEnvZR-kbSH3.*ScOmpRP.*EnvZ+kdSH3.*ScEnvZRP-kbSH3.*EnvZ+kdSH3.*EnvZD-fdE.*EnvZ+fpE;%new
        dEnvZDdt=kbSH3.*EnvZ-kdSH3.*EnvZD;%new
        dEnvZPdt=kap.*EnvZ-kad.*EnvZP-kb1.*EnvZP.*(OmpR.^2)+kd1.*EnvZPR-kbSH3.*Sc.*EnvZP+kdSH3.*ScEnvZP-kbSH3.*ScOmpR.*EnvZP+kdSH3.*ScEnvZPR;
        dEnvZPRdt=kb1.*EnvZP.*(OmpR.^2)-kd1.*EnvZPR-kpt.*EnvZPR;
        dEnvZRPdt=kpt.*EnvZPR-kd2.*EnvZRP+kb2.*EnvZ.*OmpRP-kph.*EnvZRP;
        dEnvZRdt=kph.*EnvZRP+kb3.*EnvZ.*(OmpR.^2)-kd3.*EnvZR;
        %dOmpRdt=-kb1.*EnvZP.*OmpR+kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dOmpRdt=-2.*kb1.*EnvZP.*(OmpR.^2)+2.*kd1.*EnvZPR+2.*kd3.*EnvZR-2.*kb3.*EnvZ.*(OmpR.^2)-kbLZ.*Sc.*OmpR+kdLZ.*ScOmpR-kbLZ.*ScEnvZ.*OmpR+kdLZ.*ScEnvZR-kbLZ.*ScEnvZP.*OmpR+kdLZ.*ScEnvZPR-fdO.*OmpR+fpO;
        dOmpRPdt=kd2.*EnvZRP-kb2.*EnvZ.*OmpRP-kbLZ.*Sc.*OmpRP+kdLZ.*ScOmpRP-kbLZ.*ScEnvZ.*OmpRP+kdLZ.*ScEnvZRP;
        OmpRPt=x(7)+x(15);
        dScdt=-kbSH3.*Sc.*EnvZ-kbSH3.*Sc.*EnvZP-kbLZ.*Sc.*OmpR-kbLZ.*Sc.*OmpRP+kdSH3.*ScEnvZ+kdSH3.*ScEnvZP+kdLZ.*ScOmpR+kdLZ.*ScOmpRP+kG.*kC.*OmpRPt.^2./(OmpRPt.^2+KC.^2)-kdScaf.*Sc;
        dScEnvZdt=-kap2.*ScEnvZ+kad.*ScEnvZP+kbSH3.*Sc.*EnvZ-kdSH3.*ScEnvZ-kbLZ.*ScEnvZ.*OmpR+kdLZ.*ScEnvZR-kbLZ.*ScEnvZ.*OmpRP+kdLZ.*ScEnvZRP;
        dScEnvZPdt=kap2.*ScEnvZ-kad.*ScEnvZP+kbSH3.*Sc.*EnvZP-kdSH3.*ScEnvZP-kbLZ.*ScEnvZP.*OmpR+kdLZ.*ScEnvZPR;
        dScEnvZPRdt=kbLZ.*ScEnvZP.*OmpR-kdLZ.*ScEnvZPR+kbSH3.*ScOmpR.*EnvZP-kdSH3.*ScEnvZPR-kpt2.*ScEnvZPR;%+kap2.*ScEnvZR-kad.*ScEnvZPR;%new
        dScEnvZRPdt=kbLZ.*ScEnvZ.*OmpRP-kdLZ.*ScEnvZRP+kbSH3.*ScOmpRP.*EnvZ-kdSH3.*ScEnvZRP+kpt2.*ScEnvZPR-kph.*ScEnvZRP;
        dScEnvZRdt=kbLZ.*ScEnvZ.*OmpR-kdLZ.*ScEnvZR+kbSH3.*ScOmpR.*EnvZ-kdSH3.*ScEnvZR+kph.*ScEnvZRP;%-kap2.*ScEnvZR+kad.*ScEnvZPR;%new
        dScOmpRdt=kbLZ.*Sc.*OmpR-kdLZ.*ScOmpR-kbSH3.*ScOmpR.*EnvZ+kdSH3.*ScEnvZR-kbSH3.*ScOmpR.*EnvZP+kdSH3.*ScEnvZPR;
        dScOmpRPdt=kbLZ.*Sc.*OmpRP-kdLZ.*ScOmpRP-kbSH3.*ScOmpRP.*EnvZ+kdSH3.*ScEnvZRP;
        %dGFPdt=kG.*kC.*OmpRP.^2./(OmpRP.^2+KC.^2)-kdG.*GFP;
        %dRFPdt=kR.*kF.*OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)))-kdR.*RFP;
        %;dGFPdt;dRFPdt
        EC=OmpRPt.^2./(OmpRPt.^2+KC.^2);
        EF=OmpRPt.^2./(OmpRPt.^2+KF.^2).*(1-(OmpRPt.^2./(OmpRPt.^2+KF4.^2)));
        dmGFPdt=kmGtx.*EC-kasGb.*mGFP.*asGFP+(kasGd+kDasG).*asmGFP-kDmG.*mGFP;
        dasGFPdt=kasGtx.*EF-kasGb.*mGFP.*asGFP+(kasGd+kDmG).*asmGFP-kDasG.*asGFP;
        dasmGFPdt=kasGb.*mGFP.*asGFP-(kasGd+kDmG+kDasG).*(asmGFP);
        dmRFPdt=kmRtx.*EF-kasRb.*mRFP.*asRFP+(kasRd+kDasR).*asmRFP-kDmR.*mRFP;
        dasRFPdt=kasRtx.*EC-kasRb.*mRFP.*asRFP+(kasRd+kDmR).*asmRFP-kDasR.*asRFP;
        dasmRFPdt=kasRb.*mRFP.*asRFP-(kasRd+kDmR+kDasR).*(asmRFP);
        dGFPpdt=kGtl.*mGFP-(kfG+kDG).*(GFPp);
        dGFPadt=kfG.*GFPp-kDG.*GFPa;
        dRFPpdt=kRtl.*mRFP-(kfR+kDR).*(RFPp);
        dRFPadt=kfR.*RFPp-kDR.*RFPa;
        dxdt=[dEnvZdt;dEnvZPdt;dEnvZPRdt;dEnvZRPdt;dEnvZRdt;dOmpRdt;dOmpRPdt;dScdt;dScEnvZdt;dScEnvZPdt;dScEnvZPRdt;dScEnvZRPdt;dScEnvZRdt;dScOmpRdt;dScOmpRPdt;dEnvZDdt;dmGFPdt;dasGFPdt;dasmGFPdt;dmRFPdt;dasRFPdt;dasmRFPdt;dGFPpdt;dGFPadt;dRFPpdt;dRFPadt];
    end

    function [position, isterminal, direction] = reachSS(t, X)
      Changee=abs(f(t, X));
      maxChange = max(Changee(1:15));
      %fprintf('SS \n');
      position = (maxChange < 1e-10) - 0.5;        
      isterminal = 1;        
      direction = 1;    
    end 
end

