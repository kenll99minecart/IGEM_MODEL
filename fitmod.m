function lsbeta2 = fitmod()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Testing fitting function

%save file:
x=zeros(48,2);
x1=zeros(12,2);
x1(1:3:12,2)=15;
x1(2:3:12,2)=45;
x1(3:3:12,2)=90;
x1(1:3,1)=127;
x1(4:6,1)=516;
x1(7:9,1)=639;
x1(10:12,1)=1010;
%since x1 indicates sucrose, the figure may not be accurate

x2=zeros(27,2);
n=[5,10,15,30,45,60,90,100,120];
x2(1:9,2)=n;
x2(10:18,2)=n;
x2(19:27,2)=n;
x2(1:9,1)=17e-3;
x2(10:18,1)=500e-3;
x2(19:27,1)=850e-3;
%save x2;
%30 41 58 72 82 83 85 95 90
%23 41 58 61 70 71 69 70 70
EnvZ=[12 20 18 21 23 23 23 25 30 30 41 58 72 82 83 85 95 90 5 6 7 10 11 11 11 10 10 23 41 58 61 70 71 69 70 70]';
EnvZ=((1-(EnvZ./100)).*4);
sodium = [17e-3;500e-3]; %in \mu M ;850e-3
Osm = [299e-3;1214e-3;133e-3;1728e-3];
tt=n.*60;
beta0 = 3e-7; %initial guess V1
beta1 = [1,1,1]; %initial guess V2
modelfunc = @(params, tt) constfunc(tt, params,sodium, 0);
modelfunc1 = @(params, tt)constfuncT(tt, params, sodium, 0);
modelfunc2 = @(params, tt,sodium)constfunc(tt, params, sodium, 0);
modelfuncOsm = @(params, tt)constfuncT(tt,params,Osm,0);
a=size(constfunc(tt,beta0,sodium,0));
b=size(EnvZ);
%[lsbeta2, resnorm, residual, exitflag, output, lambda, jacobian] =lsqcurvefit(modelfunc, beta0, tt, EnvZ,0,1/max(sodium));  %V1
[lsbeta2, resnorm, residual, exitflag, output, lambda, jacobian] =lsqcurvefit(modelfuncOsm, beta1, tt, EnvZ,[0,-Inf,0]); %V2
%[lsbeta2, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(modelfunc,beta0,tt,EnvZ,'Y',sodium);
%ci2 = nlparci(lsbeta2, residual, 'Jacobian', jacobian);
%beta = nlinfit(tt, EnvZ, @(params, tt) constfunc(tt, params,sodium, 0), beta0);
%[x,resnorm, ~, exitflag, output] = lsqcurvefit(modelfunc,beta0, tt, EnvZ);
%fit function
% EnvZmat=reshape(EnvZ,[9,3]);
% EnvZmat=EnvZmat';
% [fitX,fitY,fitZ]=prepareSurfaceData(tt,sodium,EnvZmat);
% [fitX,fitY];
% ft=fittype('constfuncU(x,z,params)','independent',{'x','z'},'coefficients',{'params'});
% fitmod=fit([fitX,fitY],fitZ,ft,'StartPoint',3e-9,'Lower',0);
% fitmod
tf= max(tt);
%  figure(1);
 %plot(fitmod,[fitX,fitY],fitZ);
 hold off;
 constfuncT(tf, lsbeta2, Osm, 1);  % plot the model curve with the fitted parameters
 hold on;
 time= zeros(9*4,1);
 time(1:9)=n;
 time(10:18)=n;
 time(19:27)=n;
 time(28:36)=n;
 %plot(x2(1:end-9,2).*60, EnvZ, 'kx');
 plot(time.*60, EnvZ, 'kx');% add the data points to the plot
% figure(1)
%  plot(x2(1:end-9,2).*60,constfunc(tt, x, sodium, 0));
%  hold on;
%  plot(x2(1:end-9,2).*60, EnvZ, 'kx'); 
%  hold off;
end

