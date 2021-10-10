function Ans = TCSanti2(kap,showplot,initial)
% important:  according to literature, EnvZ must be far greater than OmpR

%parameters goes here
kad = 0.0001;%0.001
kb1 = 0.5;
kd1 = 0.5;
kpt= 0.0015; %1.5/0.0015 %need to sensitivity analysis
kd2 = 0.5;
kb2 = 0.5;
kph =0.025;%0.0015e-6%wrong prev 0.05
kd3 = 0.5;
kb3 = 0.5;
tf = 100000;

KC=20e-3;
KF=1e-3;
KF4=20e-3;

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
% kG=1;
% kR=1;
% kC=0.01;
% kF=0.01;
% kdG=0.001;
% kdR=0.001;
% kdGm=1;
% kdRm=1;
% KdG=1e-2;
% KdR=1e-2;
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
%initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi];
%if(~SS)
opt=odeset('RelTol',1e-5,'NonNegative',1,'Events',@reachSS);%
[tt, NN] = ode15s(@f, [0; tf], initial,opt);

Ans=NN(end,:);
%else
%guess=ones(1,9).*7;
%ops=optimoptions(@fsolve,'FunctionTolerance',1e-20);
%ops2=optimoptions(@lsqnonlin,'FunctionTolerance',1e-20,'MaxIterations',1e10,'MaxFunctionEvaluations',1e10,'StepTolerance',1e-30,'Algorithm','trust-region-reflective');
t=1000;
%Ans=fsolve(@(x)fSS(t,x),initial,ops);
%Ans=lsqnonlin(@(x)fSS(t,x),guess,zeros(1,9),[],ops2);
%end

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
% plot(tt,NN(:,8),tt,NN(:,9),tt,NN(:,10),tt,NN(:,11));
% legend('GFPt','RFPt','GFPL','RFPL');
% xlabel('Time / s');
% ylabel('Number of moles/\muM');
figure(3);
plot(tt,NN(:,15),tt,NN(:,17));
legend('GFPL','RFPL');
xlabel('Time / s');
ylabel('Number of moles/\muM');
figure(4);
plot(tt,NN(:,7)+NN(:,6)+NN(:,5)+NN(:,4)+NN(:,3));
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
        mGFP=x(8);
        asGFP=x(9);
        asmGFP=x(10);
        mRFP=x(11);
        asRFP=x(12);
        asmRFP=x(13);
        GFPp=x(14);
        GFPa=x(15);
        RFPp=x(16);
        RFPa=x(17);
        dEnvZdt=-kap.*EnvZ+kad.*EnvZP+kd2.*EnvZRP-kb2.*EnvZ.*OmpRP+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dEnvZPdt=kap.*EnvZ-kad.*EnvZP-kb1.*EnvZP.*(OmpR.^2)+kd1.*EnvZPR;
        dEnvZPRdt=kb1.*EnvZP.*(OmpR.^2)-kd1.*EnvZPR-kpt.*EnvZPR;
        dEnvZRPdt=kpt.*EnvZPR-kd2.*EnvZRP+kb2.*EnvZ.*OmpRP-kph.*EnvZRP;
        dEnvZRdt=kph.*EnvZRP+kb3.*EnvZ.*OmpR-kd3.*EnvZR;
        dOmpRdt=-2.*kb1.*EnvZP.*(OmpR.^2)+2.*kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dOmpRPdt=kd2.*EnvZRP-kb2.*EnvZ.*OmpRP;
        %assume constant plasmid number CN
        %dGFPtdt=kC.*OmpRP.^2./(OmpRP.^2+KC.^2)-kdGm.*GFPt;
        %dRFPtdt=kF.*OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)))-kdRm.*RFPt;
        EC=OmpRP.^2./(OmpRP.^2+KC.^2);
        EF=OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)));
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
%         dGFPLdt=KdG./(KdG+RFPt).*dGFPtdt+GFPt.*(-KdG./(KdG+RFPt).^2).*dRFPtdt;
%         dRFPLdt=KdR./(KdR+GFPt).*dRFPtdt+RFPt.*(-KdR./(KdR+GFPt).^2).*dGFPtdt;
%         dGFPpdt=kG.*GFPL-kdG.*GFPp;
%         dRFPpdt=kG.*RFPL-kdR.*RFPp;
        %;dGFPdt;dRFPdt
        dxdt=[dEnvZdt;dEnvZPdt;dEnvZPRdt;dEnvZRPdt;dEnvZRdt;dOmpRdt;dOmpRPdt;dmGFPdt;dasGFPdt;dasmGFPdt;dmRFPdt;dasRFPdt;dasmRFPdt;dGFPpdt;dGFPadt;dRFPpdt;dRFPadt];
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

