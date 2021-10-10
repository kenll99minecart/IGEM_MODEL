function [] = StabilityAnalysis()
%UNTITLED Summary of this function goes here

vid=VideoWriter('scaffold.mp4','MPEG-4');
vid.FrameRate=60;
open(vid)
Vt=ones(15,15);
Ut=ones(15,15);
for i=1:max(length(OmpRP))

x=linspace(0,30,15);
y=linspace(0,30,15);
[X,Y]=meshgrid(x,y);
grid on;
TCSwscafd2
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
drawnow limitrate
F(i) = getframe(gcf);
hold off;
writeVideo(vid,F(i));
end

%close(vid);
%movie(F);


end

