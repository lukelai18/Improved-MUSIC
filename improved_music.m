p=4;%入射信号数目
M=7;%阵元
fc=1e9;%入射信号中频
DOA=[0,10,30,50];%信号入射角度
fs=3*fc;%采样频率
N=300;%采样点个数
dt=1/fs;%采样时间间隔
t=0:dt:(N-1)*dt;%时间轴
c=3e8;%波速
d=c/fc/2;%阵元间距为半波长

s1=bkskdigital([1 0 1 1 1 0 0 1],fc);
s2=qpsk(8,fc);
%s1=sqrt(2)*cos(2*pi*fc*t);%信号1
%mt=sqrt(2)*cos(2*pi*5e8*t);
%s2=mt.*cos(2*pi*fc*t);%信号2
%mt=sqrt(2)*cos(2*pi*1e7*t+pi/8);
%s3=mt.*cos(2*pi*fc*t);%信号3
%ss=[s1;s2;s3];
%s=ss(1:p,:);

%case1
ss1=[s1;s2;conj(s2);s2;conj(s2);s2;conj(s2)];
s_1=ss1(1:7,:);
DOA=[0,10,0,30,0,50,0];%信号入射角度
a1=ula(M,d,DOA(1),fc);
A(1:M,1)=a1;
A(M+1:2*M,1)=conj(a1);
for k=2:2:7
    a2=ula(M,d,DOA(k),fc);
    A(1:M,k)=a2;
    A(M+1:2*M,k)=zeros(1,M)';
    A(1:M,k+1)=zeros(1,M)';
    A(M+1:2*M,k+1)=conj(a2);
end
y1=[A*s_1;conj(A)*conj(s_1)];

%%%%%RMSE计算
k=1;                                                                    
xll=-10:1:20;                            %信噪比范围
for SNR=-0:5:30                        
    sumtheta_1=0;                        %误差
    for i=1:100                            
        xt=awgn(y1,SNR);                  
        [theta_1]=im_music(y1,SNR,length(DOA),M,N,fc,d,c); %music估计
        theta_1=sort(theta_1); 
        sumtheta_1=sumtheta_1+(((0)-theta_1(1,1))^2+(10-theta_1(1,2))^2+(30-theta_1(1,3))^2+(50-theta_1(1,4))^2);
    end      
    RMS_1(k)=sqrt(sumtheta_1/4/100);   
    k=k+1;
end


%case2
ss1=[s1;s1;s2;conj(s2);s2;conj(s2)];
s_1=ss1(1:6,:);
DOA=[0,10,30,0,50,0];%信号入射角度
for k=1:2
a1=ula(M,d,DOA(k),fc);
A(1:M,k)=a1;
A(M+1:2*M,k)=conj(a1);
end
for k=3:2:6
    a2=ula(M,d,DOA(k),fc);
    A(1:M,k)=a2;
    A(M+1:2*M,k)=zeros(1,M)';
    A(1:M,k+1)=zeros(1,M)';
    A(M+1:2*M,k+1)=conj(a2);
end
y2=[A*s_1;conj(A)*conj(s_1)];
snr=15;%信噪比
y2=awgn(y2,snr);%按信噪比SNR对数据y加相应的高斯白噪声
R=y2*y2'/N;%协方差矩阵

%case3
ss1=[s1;s1;s1;s2;conj(s2)];
s_1=ss1(1:5,:);
DOA=[0,10,30,50,0];%信号入射角度
for k=1:3
a1=ula(M,d,DOA(k),fc);
A(1:M,1)=a1;
A(M+1:2*M,1)=conj(a1);
end
    a2=ula(M,d,DOA(4),fc);
    A(1:M,4)=a2;
    A(M+1:2*M,4)=zeros(1,M)';
    A(1:M,5)=zeros(1,M)';
    A(M+1:2*M,5)=conj(a2);
y3=[A*s_1;conj(A)*conj(s_1)];
snr=15;%信噪比
y3=awgn(y3,snr);%按信噪比SNR对数据y加相应的高斯白噪声
R=y3*y3'/N;%协方差矩阵

%响应矩阵
for k=1:p
    a1=ula(M,d,DOA(k),fc);
    A(:,k)=a1;
end

for k=1:p
    a1=ula(M,d,DOA(k),fc);
    A(:,k)=a1;
end

y1=[A*s;A'*s'];
snr=15;%信噪比
y1=awgn(y1,snr);%按信噪比SNR对数据y加相应的高斯白噪声
R=y1*y1'/N;%协方差矩阵

%非圆信号估计
v=[A zeros(M,1);zeros(M,1) A'];
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

%侧向结果
for k=1:p
    [k1,k2]=max(pu);
    DOA_guiji(k)=(k2-1)-90;
    pu(k2)=0;
end
DOA_guiji