function Amount = constfuncU(tt,sodium,params)
%CONSTFUNC Summary of this function goes here
%   Detailed explanation goes here
showplot=0;
% if(min(size(sodium))>1)
% size(sodium)
% error('sodium matrix must be 1');
% end
a = params; 
%Sodium = [0.127;0.516;0.639;1.010;0.133;0.299;1.214;1.728]; % unit are in mOsm/g
%Sodium = x(:,1);
%sodium = [17e-6;500e-6;850e-6]; 
%ones(max(size(sodium)),1).*
Tfinal = max(tt);
Maxt=max(size(tt));
MaxS=max(size(sodium));
EnvZi = [4,0]; % in \mu M
kad = 0.001;
%figure(1);
%hold on;
iter=1;
Amountnull = ones(MaxS*Maxt,1);
opt=odeset('RelTol',1e-13,'MaxStep',0.5);
size(Amountnull);
for (i=1:MaxS)
[ttnew,yy] = ode45(@f, [0, 8000], EnvZi,opt);
Amountmat = interp1(ttnew,yy(:,1),tt);
Amountnull(1+Maxt*(i-1):Maxt*i)=Amountmat(:);
if(showplot)
plot(ttnew,yy(:,1));
end
iter=iter+1;
end

Amount = Amountmat;
if(showplot)
    xlabel('Time');
    ylabel('EnvZ');
end

    function dXdt = f(t, X)
        EnvZ=X(1);
        EnvZP=X(2);
        if (EnvZ<=0)
            EnvZ=0;
        end
        if(EnvZP<=0)
            EnvZP=0;
        end
        
        dEnvZdt = -a.*sodium(iter).*EnvZ+kad.*EnvZP;
        dEnvZPdt = a.*sodium(iter).*EnvZ-kad.*EnvZP;
        dXdt = [dEnvZdt;dEnvZPdt]; %V1
    end
end

