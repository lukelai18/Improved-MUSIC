%设置均匀阵列矩阵
function av=ula(m,d,theta,fc)
%m为阵元个数，d为阵元间距，theta为信号入射角度，fc为载频
%theta=30;m=3;fc=500;  %假设阵元间距为lamda/2，入射角度为30°，载频为500Hz,阵元个数为3
%设置传感器的矩阵
%b=zeros(m,2);
%for n=1:m
%b(n,1)=(n-1)*100;
%b(n,2)=0;
%end
lamda=3e8/fc;%波长
%av=zeros(m,1);
a=exp(-1i*2*pi*d*sind(theta)/lamda);
for n=1:m
    %av(n)=exp(j*2*pi*(b(n,1)*sind(theta)+b(n,2)*cosd(theta))/lamda);
    av(n)=a.^(n-1);
end
av=av';
