p=5;%�����ź���Ŀ
M=7;%��Ԫ
fc=1e9;%�����ź���Ƶ
DOA=[-20,-10,0,20,40];%DOA=[-30,-10,0,30,60];DOA=[-15,-7,0,15,30];DOA=[-10,-5,0,10,20];DOA=[-20,-10,0,20,40],DOA=[-25,-10,0,25,50];
%�ı��ź�����Ƕȼ�࣬��̽��MUSIC�㷨�������ź�����Ƕȼ��Ĺ�ϵ
fs=3*fc;%����Ƶ��
N=512;%������������������������������£��ֱ���Ϊ64��256��512��������ʵ��
dt=1/fs;%����ʱ����
t=0:dt:(N-1)*dt;%ʱ����
c=3e8;%����
d=c/fc/2;%��Ԫ���Ϊ�벨��

s1=sqrt(2)*cos(2*pi*fc*t);%�ź�1
mt=sqrt(2)*cos(2*pi*5e8*t);
s2=mt.*cos(2*pi*fc*t);%�ź�2
mt=sqrt(2)*cos(2*pi*1e7*t+pi/8);
s3=mt.*cos(2*pi*fc*t);%�ź�3
mt=sqrt(2)*cos(2*pi*2*fc*t+pi/4);
s4=mt.*cos(2*pi*fc*t);%�ź�4
mt=sqrt(2)*cos(2*pi*2*fc*t);
s5=mt.*cos(2*pi*fc*t);%�ź�5

ss=[s1;s2;s3;s4;s5];
s=ss(1:p,:);

%��Ӧ����
for k=1:p
    a1=ula(M,d,DOA(k),fc);
    A(:,k)=a1;
end

y=A*s;
snr=15;%����ȣ��������������������£��ֱ���Ϊ0��5��15��������ʵ��
y=awgn(y,snr);%�������SNR������y����Ӧ�ĸ�˹������
R=y*y'/N;%Э�������

%�����ӿռ����
pg=p;
[v,dd]=eig(R);
if(dd(1,1)>dd(2,2))
    Un=v(:,pg+1:M);
else
    Un=v(:,1:(M-pg));
end

%����music��
do=-90:90;
pu=zeros(1,length(do));
kg=1;
for k=-90:90
    a=zeros(M,1);
    for kk=1:M
        a(kk,1)=exp(-j*2*pi*fc*(kk-1)*d*sin(k/180*pi)/c);
    end
    pu(1,kg)=1/(a'*Un*Un'*a);%����
    kg=kg+1;
end
plot(do,10*log10(abs(pu)),'r-','linewidth',2);
grid on;
title('MUSIC����');
xlabel('���﷽��');
ylabel('MUSIC��');

%������
for k=1:p
    [k1,k2]=max(pu);
    DOA_guiji(k)=(k2-1)-90;
    pu(k2)=0;
end
DOA_guiji