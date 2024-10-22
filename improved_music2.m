p=6;%入射信号数目
M=6;%阵元
fc=1e9;%入射信号中频
%DOA=[-40,-20,0,10,30,50];%信号入射角度
fs=5*fc;%采样频率
N=5000;%采样点个数
dt=1/fs;%采样时间间隔
t=0:dt:(N-1)*dt;%时间轴
c=3e8;%波速
d=c/fc/2;%阵元间距为半波长

s1=bkskdigital(round(rand(1,8)),fc); %输入信号也可能会有问题
s2=bkskdigital(round(rand(1,8)),fc);
s3=bkskdigital(round(rand(1,8)),fc);
s4=bkskdigital(round(rand(1,8)),fc);
s5=qpsk(8,fc); 
s6=qpsk(8,fc); 

%case1
ss1=[s1;s2;s3;s4;s5;conj(s5);s6;conj(s6)];  %非圆度信号与圆度信号都是一致的，可能需要修改
s_1=ss1(1:8,:);
DOA=[-40,-20,0,10,30,0,50,0];%信号入射角度

v1=zeros(1,8);
for k=1:8
    v1(k)=1/cosd(DOA(k));
end
b1=diag(v1);%文献中第一个B的构造（信道对信号的影响）
v2=zeros(1,8);
for k=1:4
    v2(k)=1/cosd(DOA(k));
end
for k=5:8
    v2(k)=3; %后面可以考虑换成1，调试中发现后面的值均变成了3
end
b2=diag(v2);%文献中第二个B的构造（表示相位旋转，对圆度信号不成立）
b3=b1*b2;%表示对信号的总影响

for k=1:4
a1=ula(M,d,DOA(k),fc);
A(1:M,k)=a1*b3(k,k);
A(M+1:2*M,k)=conj(a1)*conj(b3(k,k));
end
for k=5:2:8
    a2=ula(M,d,DOA(k),fc);
    A(1:M,k)=a2*b3(k,k);
    A(M+1:2*M,k)=zeros(1,M)';
    A(1:M,k+1)=zeros(1,M)';
    A(M+1:2*M,k+1)=conj(a2)*conj(b3(k,k));
end
y1=A*s_1;
snr=20;%信噪比
y1=awgn(y1,snr);%按信噪比SNR对数据y加相应的高斯白噪声
R=y1*y1'/N;%协方差矩阵

%噪声子空间估计
pg=length(DOA);
[v,dd]=eig(R);
if(dd(1,1)>dd(2,2))
    Un=v(:,pg+1:2*M);
else
    Un=v(:,1:(2*M-pg));
end

%一半的噪声子空间
for k=1:2
    Un1(:,k)=Un(:,k);
end

%计算music谱
do=-90:90;
pu=zeros(1,length(do));
kg=1;
kt=2;
b=zeros(M,2);
%for n=1:M
%b(n,1)=(n-1)*100;
%b(n,2)=0;
%end
for k=-90:90
    a=zeros(2*M,2);
    for kk=1:M
        %a(kk,1)=exp(1i*2*pi*(b(kk,1)*sind(k)+b(kk,2)*cosd(k))/d);
        a(kk,1)=exp(-1i*2*pi*fc*(kk-1)*d*sind(-k)/c); %备注：与原模型相比。把sind里的k改成了负数
        a(kk,2)=0;
    end
     for kk=M+1:2*M
       a(kk,1)=0;
       a(kk,2)=conj(exp(-1i*2*pi*fc*(kk-M-1)*d*sind(-k)/c));
       %a(kk,2)=conj(exp(1i*2*pi*(b(kk-M,1)*sind(k)+b(kk-M,2)*cosd(k))/d));
     end
     e=zeros(2*M,1);
     for kk=1:M
        e(kk,1)=0;
    end
     for kk=M+1:2*M
       e(kk,1)=conj(exp(-1i*2*pi*fc*(kk-M-1)*d*sind(-k)/c));%备注：与原模型相比。把sind里的k改成了负数
     end
    pu(1,kg)=1/det(a'*Un*Un'*a);%算法[29]谱线
    pu1(1,kg)=1/det(e'*Un1*Un1'*e);%算法[35]谱线
    kg=kg+1;
end
plot(do,10*log10(abs(pu)),'r-','linewidth',2);hold on;
plot(do,10*log10(abs(pu1)),'b-','linewidth',2)
grid on;
title('MUSIC测向');
xlabel('波达方向');
ylabel('MUSIC谱');
