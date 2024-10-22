p=6;%�����ź���Ŀ
M=6;%��Ԫ
fc=1e9;%�����ź���Ƶ
%DOA=[-40,-20,0,10,30,50];%�ź�����Ƕ�
fs=5*fc;%����Ƶ��
N=5000;%���������
dt=1/fs;%����ʱ����
t=0:dt:(N-1)*dt;%ʱ����
c=3e8;%����
d=c/fc/2;%��Ԫ���Ϊ�벨��

s1=bkskdigital(round(rand(1,8)),fc); %�����ź�Ҳ���ܻ�������
s2=bkskdigital(round(rand(1,8)),fc);
s3=bkskdigital(round(rand(1,8)),fc);
s4=bkskdigital(round(rand(1,8)),fc);
s5=qpsk(8,fc); 
s6=qpsk(8,fc); 

%case1
ss1=[s1;s2;s3;s4;s5;conj(s5);s6;conj(s6)];  %��Բ���ź���Բ���źŶ���һ�µģ�������Ҫ�޸�
s_1=ss1(1:8,:);
DOA=[-40,-20,0,10,30,0,50,0];%�ź�����Ƕ�

v1=zeros(1,8);
for k=1:8
    v1(k)=1/cosd(DOA(k));
end
b1=diag(v1);%�����е�һ��B�Ĺ��죨�ŵ����źŵ�Ӱ�죩
v2=zeros(1,8);
for k=1:4
    v2(k)=1/cosd(DOA(k));
end
for k=5:8
    v2(k)=3; %������Կ��ǻ���1�������з��ֺ����ֵ�������3
end
b2=diag(v2);%�����еڶ���B�Ĺ��죨��ʾ��λ��ת����Բ���źŲ�������
b3=b1*b2;%��ʾ���źŵ���Ӱ��

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
snr=20;%�����
y1=awgn(y1,snr);%�������SNR������y����Ӧ�ĸ�˹������
R=y1*y1'/N;%Э�������

%�����ӿռ����
pg=length(DOA);
[v,dd]=eig(R);
if(dd(1,1)>dd(2,2))
    Un=v(:,pg+1:2*M);
else
    Un=v(:,1:(2*M-pg));
end

%һ��������ӿռ�
for k=1:2
    Un1(:,k)=Un(:,k);
end

%����music��
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
        a(kk,1)=exp(-1i*2*pi*fc*(kk-1)*d*sind(-k)/c); %��ע����ԭģ����ȡ���sind���k�ĳ��˸���
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
       e(kk,1)=conj(exp(-1i*2*pi*fc*(kk-M-1)*d*sind(-k)/c));%��ע����ԭģ����ȡ���sind���k�ĳ��˸���
     end
    pu(1,kg)=1/det(a'*Un*Un'*a);%�㷨[29]����
    pu1(1,kg)=1/det(e'*Un1*Un1'*e);%�㷨[35]����
    kg=kg+1;
end
plot(do,10*log10(abs(pu)),'r-','linewidth',2);hold on;
plot(do,10*log10(abs(pu1)),'b-','linewidth',2)
grid on;
title('MUSIC����');
xlabel('���﷽��');
ylabel('MUSIC��');
