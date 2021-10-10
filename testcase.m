function x = testcase()
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
%initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi;GFPi;RFPi];
       ops2=optimoptions(@lsqnonlin,'FunctionTolerance',1e-10,'MaxIterations',1e10,'MaxFunctionEvaluations',1e10,'StepTolerance',1e-30);
       [h,resnorm,residual,exitflag,output,lambda,jacobian] =lsqnonlin(@(x)f(0,x),ones(1,26).*1e2,zeros(1,26));%ones(1,17).*1
       [k,resnorm2,residual2,exitflag2,output2,lambda2,jacobian2] =lsqnonlin(@(x)f(0,x),ones(1,26).*1,zeros(1,26));
       
       kaprangeSS=linspace(0,0.01,101);
       RFPrange=linspace(0,100,101);
       GFPrange=linspace(0,100,101);
       [X,Y]=meshgrid(RFPrange,GFPrange);
       RFPpt=h(25);
       R=kfR.*RFPpt-kDR.*X;
       G=
       quiver(X,Y,
       x=linspace(0,1000,50);
       y=linspace(0,1000,50);
       %[x,y]=meshgrid(0:10:101,0:10:101);
       h=h';
       k=k';
       [U,V]=meshgrid(x,y);
       U
       V
       h
       k
       size(V)
       size(h)
       t=0;
       Gcood=zeros(50,50);
       Rcood=zeros(50,50);
       Gvec=zeros(50,50);
       Rvec=zeros(50,50);
       
       for i=1:50
           for j=1:50
       w=([h(24),k(24);h(26),k(26)])\[x(i);y(j)];
           omega=h*w(1)+k*w(2);
           omega
           Gcood(i,j)=x(i);
           Rcood(i,j)=y(j);
           mat=f(t,omega);
           Gvec(i,j)=mat(24);
           Rvec(i,j)=mat(26);
           end
       end
       quiver(Gcood,Rcood,Gvec,Rvec);
       h
       k
       x;
       [V,ansf]=eig(full(jacobian));
       ansf
       opt=odeset('RelTol',1e-5,'NonNegative',1,'Events',@reachSS);%
     tf=1000000;
     resultsOmpRP=ones(1,200);
     initialv=linspace(0,2e7,200);
%      for i=1:200
%      initial(6)=initialv(i);
%      [tt, NN] = ode15s(@f, [0; tf], initial,opt);
%     resultsOmpRP(i)=NN(end,7);
%      end
%      plot(initialv,resultsOmpRP);
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
        dEnvZRdt=kph.*EnvZRP+kb3.*EnvZ.*OmpR-kd3.*EnvZR;
        %dOmpRdt=-kb1.*EnvZP.*OmpR+kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dOmpRdt=-2.*kb1.*EnvZP.*(OmpR.^2)+2.*kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR-kbLZ.*Sc.*OmpR+kdLZ.*ScOmpR-kbLZ.*ScEnvZ.*OmpR+kdLZ.*ScEnvZR-kbLZ.*ScEnvZP.*OmpR+kdLZ.*ScEnvZPR-fdO.*OmpR+fpO;
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
      maxChange = max(Changee(1:7));
      %fprintf('SS \n');
      position = (maxChange < 1e-10) - 0.5;        
      isterminal = 1;        
      direction = 1;    
    end
end