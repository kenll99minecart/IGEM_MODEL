function ans = OmpR(minR,maxR)
%SS for OmpR versus 

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
KdG=1e0;
KdR=1e0;
OmpRP=linspace(minR,maxR,100);
% ansNorm=zeros(2,100);
% iteration=1;
% ops=optimoptions(@fsolve,'FunctionTolerance',1e-20);
% for i = 1:100
% ansNorm(:,2)=fsolve(@(G,R) fnorm(G,R),100,100,ops);
% iteration=iteration+1;
% end

%ansAnti=fsovle(@(x) fanti(x),initial);

GFP=kG.*kC.*OmpRP.^2./(OmpRP.^2+KC.^2)./(kdG.*kdGm);
RFP=kR.*kF.*OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)))./(kdR.*kdRm);

GFPm=kC.*OmpRP.^2./(OmpRP.^2+KC.^2)./kdG;
RFPm=kF.*OmpRP.^2./(OmpRP.^2+KF.^2).*(1-(OmpRP.^2./(OmpRP.^2+KF4.^2)))./kdR;
antiG=RFPm;
antiR=GFPm;
leftGFP=GFPm.*KdG./(KdG+antiG);
leftRFP=RFPm.*KdR./(KdR+antiR);

GFPanti=kG.*leftGFP./kdGm;
RFPanti=kR.*leftRFP./kdRm;
%assume fast binding;
%0=kC.*OmpRP.^2./(OmpRP.^2+KC.^2)-kdGm.*GFPm-kb.*GFPm.*antiG+kf.*;
ans=fzero(@(OmpRPi) kF.*OmpRPi.^2./(OmpRPi.^2+KF.^2).*(1-(OmpRPi.^2./(OmpRPi.^2+KF4.^2)))./kdRm-kC.*OmpRPi.^2./(OmpRPi.^2+KC.^2)./kdGm,0.02)

figure(1);
xlabel('OmpRP / \muM');
ylabel('FPs / \muM');
plot(OmpRP,GFP,'g-',OmpRP,RFP,'r-',OmpRP,GFPanti,'b-',OmpRP,RFPanti,'c-');
legend('GFP','RFP','GFP with antisense','RFP with antisense')

figure(2);
plot(RFP,GFP,'Color',[0 1 0]);
xlabel('RFP');
ylabel('GFP')


figure(3);
plot(RFPanti,GFPanti,'b-');
xlabel('RFP');
ylabel('GFP')

%F(max(length(OmpRP))) = struct('cdata',[],'colormap',[]);
cap=figure(4);

%vid=VideoWriter('v2.mp4','MPEG-4');
%vid.FrameRate=60;
%open(vid)
Vt=ones(15,15);
Ut=ones(15,15);
for i=1:max(length(OmpRP))
x=linspace(0,30,15);
y=linspace(0,30,15);
[X,Y]=meshgrid(x,y);
grid on;
V=kG.*kC.*OmpRP(i).^2./(OmpRP(i).^2+KC.^2)-kdG.*Y;
U=kR.*kF.*OmpRP(i).^2./(OmpRP(i).^2+KF.^2).*(1-(OmpRP(i).^2./(OmpRP(i).^2+KF4.^2)))-kdR.*X;
Vt=Vt+V;
Ut=Ut+U;
GFPi=kG.*kC.*OmpRP(i).^2./(OmpRP(i).^2+KC.^2)./(kdG.*kdGm);
RFPi=kR.*kF.*OmpRP(i).^2./(OmpRP(i).^2+KF.^2).*(1-(OmpRP(i).^2./(OmpRP(i).^2+KF4.^2)))./(kdR.*kdRm);

quiver(X,Y,U,V,'r','LineWidth',0.9);
axis equal
hold on;
plot(RFP,GFP,'b-',RFPi,GFPi,'ro');
xlabel('RFP');
ylabel('GFP');
%drawnow limitrate
%F(i) = getframe(gcf);
hold off;
%writeVideo(vid,F(i));
end

%close(vid);
%movie(F);

quiver(X,Y,Ut,Vt,'m');
hold on;
plot(RFP,GFP,'b-');
xlabel('RFP');
ylabel('GFP');

figure(5);
 i=length(OmpRP)/2;
 Vc=kG.*kC.*OmpRP(i).^2./(OmpRP(i).^2+KC.^2)-kdG.*Y;
 Uc=kR.*kF.*OmpRP(i).^2./(OmpRP(i).^2+KF.^2).*(1-(OmpRP(i).^2./(OmpRP(i).^2+KF4.^2)))-kdR.*X;
 quiver(X,Y,Uc,zeros(15,15))
 hold on;
 %plot()
 %quiver(X,Y,zeros(15,15),Vc);
% x=linspace(0,30,15);
% y=linspace(0,30,15);
% [Xn,Yn]=meshgrid(x,y);
% OmpRPs=zeros(1,15);
% for i=1:15
% %ops=optimset(@fzero,'FunctionTolerance',1e-20);
% OmpRPs(i)=fzero(@(OmpRPi)kG.*kC.*OmpRPi.^2./(OmpRPi.^2+KC.^2)-kdG.*y(i),0.00005);
% end
% RFPf=kR.*kF.*OmpRPs.^2./(OmpRPs.^2+KF.^2).*(1-(OmpRPs.^2./(OmpRPs.^2+KF4.^2)))./kdR;
% plot(y,RFPf,'b-');
% 
% OmpRPf=0;
% opt=odeset('RelTol',1e-7,'NonNegative',1,'Events',@reachSS);
% tf=3000;
% OmpShort=linspace(minR,maxR,100);
% resOmp=zeros(4,100);
% for i=1:100
% OmpRPf=OmpShort(i);
% [tt, NN] = ode45(@fv, [0; tf],[0,0,0,0],opt);
% resOmp(:,i)=NN(end,:);
% end 

% figure(6);
% xlabel('OmpRP / \muM');
% ylabel('FPs / \muM');
% plot(OmpRP,GFP,'g-',OmpRP,RFP,'r-',OmpRP,GFPanti,'b-',OmpRP,RFPanti,'c-',OmpShort,resOmp(3,:),OmpShort,resOmp(4,:));
% legend('GFP','RFP','GFP with antisense','RFP with antisense','GFP with antisenseODE','RFP with antisenseODE')
% 
% figure(7);
% xlabel('OmpRP / \muM');
% ylabel('FPs / \muM');
% plot(OmpShort,resOmp(3,:),OmpShort,resOmp(4,:));
% legend('GFPL','RFPL');

OmpRPf=OmpShort(48);
[tt, NN] = ode45(@fv, [0; tf],[0,0,0,0],opt);
figure(8);
plot(tt,NN(:,3),tt,NN(:,4));

figure(9);
xlabel('OmpRP / \muM');
ylabel('FPs / \muM');
plot(OmpRP,GFP,'g-',OmpRP,RFP,'r-');
legend('GFP','RFP')

    function dxdt=fv(t,x)
        GFPt=x(1);
        RFPt=x(2);
        GFPL=x(3);
        RFPL=x(4);
        dGFPtdt=kG.*kC.*OmpRPf.^2./(OmpRPf.^2+KC.^2)-kdG.*GFPt;
        dRFPtdt=kR.*kF.*OmpRPf.^2./(OmpRPf.^2+KF.^2).*(1-(OmpRPf.^2./(OmpRPf.^2+KF4.^2)))-kdR.*RFPt;
        dGFPLdt=KdG./(KdG+RFPt).*dGFPtdt+GFPt.*(-KdG./(KdG+RFPt).^2).*dRFPtdt;
        dRFPLdt=KdR./(KdR+GFPt).*dRFPtdt+RFPt.*(-KdR./(KdR+GFPt).^2).*dGFPtdt;
    dxdt=[dGFPtdt;dRFPtdt;dGFPLdt;dRFPLdt];
    end  
    function [position, isterminal, direction] = reachSS(t, X)
      Changee=abs(fv(t, X));
      maxChange = max(Changee(1:4));    
      position = (maxChange < 1e-5) - 0.5;        
      isterminal = 1;        
      direction = 1;    
    end  
end

