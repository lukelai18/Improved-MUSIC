p=5;%入射信号数目
M=7;%阵元
fc=1e9;%入射信号中频
DOA=[-20,-10,0,20,40];%DOA=[-30,-10,0,30,60];DOA=[-15,-7,0,15,30];DOA=[-10,-5,0,10,20];DOA=[-20,-10,0,20,40],DOA=[-25,-10,0,25,50];
%改变信号入射角度间距，以探究MUSIC算法精度与信号入射角度间距的关系
fs=3*fc;%采样频率
N=512;%采样点个数，在其他条件不变的情况下，分别设为64、256和512，并进行实验
dt=1/fs;%采样时间间隔
t=0:dt:(N-1)*dt;%时间轴
c=3e8;%波速
d=c/fc/2;%阵元间距为半波长

s1=sqrt(2)*cos(2*pi*fc*t);%信号1
mt=sqrt(2)*cos(2*pi*5e8*t);
s2=mt.*cos(2*pi*fc*t);%信号2
mt=sqrt(2)*cos(2*pi*1e7*t+pi/8);
s3=mt.*cos(2*pi*fc*t);%信号3
mt=sqrt(2)*cos(2*pi*2*fc*t+pi/4);
s4=mt.*cos(2*pi*fc*t);%信号4
mt=sqrt(2)*cos(2*pi*2*fc*t);
s5=mt.*cos(2*pi*fc*t);%信号5

ss=[s1;s2;s3;s4;s5];
s=ss(1:p,:);

%响应矩阵
for k=1:p
    a1=ula(M,d,DOA(k),fc);
    A(:,k)=a1;
end

y=A*s;
snr=15;%信噪比，在其他条件不变的情况下，分别设为0、5和15，并进行实验
y=awgn(y,snr);%按信噪比SNR对数据y加相应的高斯白噪声
R=y*y'/N;%协方差矩阵

%噪声子空间估计
pg=p;
[v,dd]=eig(R);
if(dd(1,1)>dd(2,2))
    Un=v(:,pg+1:M);
else
    Un=v(:,1:(M-pg));
end

%计算music谱
do=-90:90;
pu=zeros(1,length(do));
kg=1;
for k=-90:90
    a=zeros(M,1);
    for kk=1:M
        a(kk,1)=exp(-j*2*pi*fc*(kk-1)*d*sin(k/180*pi)/c);
    end
    pu(1,kg)=1/(a'*Un*Un'*a);%谱线
    kg=kg+1;
end
plot(do,10*log10(abs(pu)),'r-','linewidth',2);
grid on;
title('MUSIC测向');
xlabel('波达方向');
ylabel('MUSIC谱');

%测向结果
for k=1:p
    [k1,k2]=max(pu);
    DOA_guiji(k)=(k2-1)-90;
    pu(k2)=0;
end
DOA_guiji