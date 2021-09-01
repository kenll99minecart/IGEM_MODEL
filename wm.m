function Ans = wm(k)
%param list

tf = 100;
NaIntra=0.1;
NaExtra=0.5;
Na0=[NaIntra;NaExtra];
opt=odeset('RelTol',1e-5,'Events',@reachSS);
[tt, NN] = ode45(@f, [0; tf], Na0,opt);
Ans=tt(end,1);
hold on;
figure(1);
plot(tt,NN(:,1),'b-');
plot(tt,NN(:,2),'r-');
xlabel('Time / s');
ylabel('Number of moles/mol');
legend('Intracellur Na','Extracellular Na');
grid on;
hold off;

    function dxdt=f(t,x)
        NaI=x(1);
        NaE=x(2);
        if(NaE<NaI)
        dNaIdt=0;
        dNaEdt=0;
        else
        dNaIdt=k.*abs((NaE-NaI)).*t;
        dNaEdt=-k.*abs((NaE-NaI)).*t; 
        end
        %rate of change of susceptible cell
        dxdt=[dNaIdt;dNaEdt];
    end

    function [position, isterminal, direction] = reachSS(t, X)
      maxChange = max(abs(f(t, X)));      
      position = (maxChange < 1e-2) - 0.5;        
      isterminal = 1;        
      direction = 1;    
    end  
end

