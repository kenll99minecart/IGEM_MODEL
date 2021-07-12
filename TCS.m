function Ans = TCS(kap)
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
tf = 10000;

%EnvZPR= EnvZP.OmpR etc.
EnvZPi=0;
EnvZPRi=0;
EnvZRPi=0;
EnvZRi=0;
OmpRPi=0;
EnvZi=0.17;
OmpRi=6;
initial = [EnvZi;EnvZPi;EnvZPRi;EnvZRPi;EnvZRi;OmpRi;OmpRPi];
opt=odeset('RelTol',1e-10,'Events',@reachSS,'NonNegative',1);
[tt, NN] = ode45(@f, [0; tf], initial,opt);

Ans=tt(end,:);
hold on;
figure(1);
plot(tt,NN(:,1));
plot(tt,NN(:,2));
plot(tt,NN(:,3));
plot(tt,NN(:,4));
plot(tt,NN(:,5));
plot(tt,NN(:,6));
plot(tt,NN(:,7));
xlabel('Time / s');
ylabel('Number of moles/\mumol');
legend('EnvZ','EnvZP','EnvZP.OmpRP','EnvZ.OmpRP','EnvZ.OmpR','OmpR','OmpRP');
grid on;
hold off;

    function dxdt=f(t,x)
        EnvZ=x(1);
        EnvZP=x(2);
        EnvZPR=x(3);
        EnvZRP=x(4);
        EnvZR=x(5);
        OmpR=x(6);
        OmpRP=x(7);
        dEnvZdt=-kap.*EnvZ+kad.*EnvZP+kd2.*EnvZRP-kb2.*EnvZ.*OmpRP+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dEnvZPdt=kap.*EnvZ-kad.*EnvZP-kb1.*EnvZP.*OmpR+kd1.*EnvZPR;
        dEnvZPRdt=kb1.*EnvZP.*OmpR-kd1.*EnvZPR-kpt.*EnvZPR;
        dEnvZRPdt=kpt.*EnvZPR-kd2.*EnvZRP+kb2.*EnvZ.*OmpRP-kph.*EnvZRP;
        dEnvZRdt=kph.*EnvZRP+kb3.*EnvZ.*OmpR-kd3.*EnvZR;
        dOmpRdt=-kb1.*EnvZP.*OmpR+kd1.*EnvZPR+kd3.*EnvZR-kb3.*EnvZ.*OmpR;
        dOmpRPdt=kd2.*EnvZRP-kb2.*EnvZ.*OmpRP;
        dxdt=[dEnvZdt;dEnvZPdt;dEnvZPRdt;dEnvZRPdt;dEnvZRdt;dOmpRdt;dOmpRPdt];
    end

    function [position, isterminal, direction] = reachSS(t, X)
      maxChange = max(abs(f(t, X)));      
      position = (maxChange < 1e-5) - 0.5;        
      isterminal = 1;        
      direction = 1;    
    end  
end

