function Ans = TCSwscar(kap,showplot,initial)
% important:  according to literature, EnvZ must be far greater than OmpR

%parameters goes here
kad = 0.001;
kb1 = 0.5;
kd1 = 0.5;
kpt= 1.5;
kd2 = 0.5;
kb2 = 0.5;
kph =0.05;
kd3 = 0.5;
kb3 = 0.5;
tf = 2000;

kap2=kap;
kpt2=1021; % ignore backwards cognate HK-RR phosphorylation (no phosphate transfer backwards)
kbLZ=0.1;
kdLZ=0.1;
kbSH3=0.1;
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
% initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi;Sci;ScEnvZPi;ScEnvZPRi;ScEnvZRPi;ScEnvZRi;ScOmpRPi;ScEnvZi;ScOmpRi];
%'RelTol',1e-10,
opt=odeset('Events',@reachSS,'NonNegative',1);
[tt, NN] = ode45(@f,[0; tf], initial,opt);

Ans=NN(end,:);
if(showplot)
% hold on;
% figure(1);
% plot(tt,NN(:,1));
% plot(tt,NN(:,2));
% plot(tt,NN(:,3));
% plot(tt,NN(:,4));
% plot(tt,NN(:,5));
% plot(tt,NN(:,6));
% plot(tt,NN(:,7));
% plot(tt,NN(:,8));
% plot(tt,NN(:,9));
% plot(tt,NN(:,10));
% plot(tt,NN(:,11));
% plot(tt,NN(:,12));
% plot(tt,NN(:,13));
% plot(tt,NN(:,14));
% plot(tt,NN(:,15));
% xlabel('Time / s');
% ylabel('Number of moles/\muM');
%  legend('EnvZ','EnvZP','EnvZP.OmpRP','EnvZ.OmpRP','EnvZ.OmpR','OmpR','OmpRP','Sc','ScEnvZ','ScEnvZP','ScEnvZPR','ScEnvZRP','ScEnvZR','ScOmpR','ScOmpRP');
% grid on;
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
        %TCS system
        dEnvZdt=-kap.*EnvZ+kad.*EnvZP+kd2.*EnvZRP-kb2.*EnvZ.*OmpRP+kd3.*EnvZR-kb3.*EnvZ.*OmpR-kbSH3.*Sc.*EnvZ+kdSH3.*ScEnvZ-kbSH3.*ScOmpR.*EnvZ+kdSH3.*ScEnvZR-kbSH3.*ScOmpRP.*EnvZ+kdSH3.*ScEnvZRP;
        dEnvZPdt=kap.*EnvZ-kad.*EnvZP-kb1.*EnvZP.*OmpR+kd1.*EnvZPR-kbSH3.*Sc.*EnvZP+kdSH3.*ScEnvZP-kbSH3.*ScOmpR.*EnvZP+kdSH3.*ScEnvZPR;
        dEnvZPRdt=kb1.*EnvZP.*OmpR-kd1.*EnvZPR-kpt.*EnvZPR;
        dEnvZRPdt=kpt.*EnvZPR-kd2.*EnvZRP+kb2.*EnvZ.*OmpRP-kph.*EnvZRP;
        dEnvZRdt=kph.*EnvZRP+kb3.*EnvZ.*OmpR-kd3.*EnvZR;
        dOmpRdt=-kb1.*EnvZP.*OmpR+kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR-kbLZ.*Sc.*OmpR+kdLZ.*ScOmpR-kbLZ.*ScEnvZ.*OmpR+kdLZ.*ScEnvZR-kbLZ.*ScEnvZP.*OmpR+kdLZ.*ScEnvZPR;
        dOmpRPdt=kd2.*EnvZRP-kb2.*EnvZ.*OmpRP-kbLZ.*Sc.*OmpRP+kdLZ.*ScOmpRP-kbLZ.*ScEnvZ.*OmpRP+kdLZ.*ScEnvZRP;
        dScdt=-kbSH3.*Sc.*EnvZ-kbSH3.*Sc.*EnvZP-kbLZ.*Sc.*OmpR-kbLZ.*Sc.*OmpRP+kdSH3.*ScEnvZ+kdSH3.*ScEnvZP+kdLZ.*ScOmpR+kdLZ.*ScOmpRP;
        dScEnvZdt=-kap2.*ScEnvZ+kad.*ScEnvZP+kbSH3.*Sc.*EnvZ-kdSH3.*ScEnvZ-kbLZ.*ScEnvZ.*OmpR+kdLZ.*ScEnvZR-kbLZ.*ScEnvZ.*OmpRP+kdLZ.*ScEnvZRP;
        dScEnvZPdt=kap2.*ScEnvZ-kad.*ScEnvZP+kbSH3.*Sc.*EnvZP-kdSH3.*ScEnvZP-kbLZ.*ScEnvZP.*OmpR+kdLZ.*ScEnvZPR;
        dScEnvZPRdt=kbLZ.*ScEnvZP.*OmpR-kdLZ.*ScEnvZPR+kbSH3.*ScOmpR.*EnvZP-kdSH3.*ScEnvZPR-kpt2.*ScEnvZPR;
        dScEnvZRPdt=kbLZ.*ScEnvZ.*OmpRP-kdLZ.*ScEnvZRP+kbSH3.*ScOmpRP.*EnvZ-kdSH3.*ScEnvZRP+kpt2.*ScEnvZPR-kph.*ScEnvZRP;
        dScEnvZRdt=kbLZ.*ScEnvZ.*OmpR-kdLZ.*ScEnvZR+kbSH3.*ScOmpR.*EnvZ-kdSH3.*ScEnvZR+kph.*ScEnvZRP;
        dScOmpRdt=kbLZ.*Sc.*OmpR-kdLZ.*ScOmpR-kbSH3.*ScOmpR.*EnvZ+kdSH3.*ScEnvZR-kbSH3.*ScOmpR.*EnvZP+kdSH3.*ScEnvZPR;
        dScOmpRPdt=kbLZ.*Sc.*OmpRP-kdLZ.*ScOmpRP-kbSH3.*ScOmpRP.*EnvZ+kdSH3.*ScEnvZRP;
        dxdt=[dEnvZdt;dEnvZPdt;dEnvZPRdt;dEnvZRPdt;dEnvZRdt;dOmpRdt;dOmpRPdt;dScdt;dScEnvZdt;dScEnvZPdt;dScEnvZPRdt;dScEnvZRPdt;dScEnvZRdt;dScOmpRdt;dScOmpRPdt];
    end

    function [position, isterminal, direction] = reachSS(t, X)
      maxChange = max(abs(f(t, X)));      
      position = (maxChange < 1e-3) - 0.5;        
      isterminal = 1;        
      direction = 1;    
    end  
end

