function [] = SAfinal()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
table=Morris();
data=table(:,2:end);
data_num=table(:,1);
data(data==-1)=0;
mean_num=zeros(10,1);
sd_num=zeros(10,1);
for i=1:size(data,1)
mean_num(i)=mean(abs(nonzeros(data(i,:))));
sd_num(i)=std(nonzeros(data(i,:)));
end
figure(2);
plot(mean_num,sd_num,'x');
text(mean_num,sd_num,num2str(data_num));
xlabel('\mu*');
ylabel('\sigma');
end

