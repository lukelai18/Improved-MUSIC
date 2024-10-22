p=4;%�����ź���Ŀ
M=7;%��Ԫ
fc=1e9;%�����ź���Ƶ
DOA=[0,10,30,50];%�ź�����Ƕ�
fs=3*fc;%����Ƶ��
N=300;%���������
dt=1/fs;%����ʱ����
t=0:dt:(N-1)*dt;%ʱ����
c=3e8;%����
d=c/fc/2;%��Ԫ���Ϊ�벨��

s1=bkskdigital([1 0 1 1 1 0 0 1],fc);
s2=qpsk(8,fc);
%s1=sqrt(2)*cos(2*pi*fc*t);%�ź�1
%mt=sqrt(2)*cos(2*pi*5e8*t);
%s2=mt.*cos(2*pi*fc*t);%�ź�2
%mt=sqrt(2)*cos(2*pi*1e7*t+pi/8);
%s3=mt.*cos(2*pi*fc*t);%�ź�3
%ss=[s1;s2;s3];
%s=ss(1:p,:);

%case1
ss1=[s1;s2;conj(s2);s2;conj(s2);s2;conj(s2)];
s_1=ss1(1:7,:);
DOA=[0,10,0,30,0,50,0];%�ź�����Ƕ�
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

%%%%%RMSE����
k=1;                                                                    
xll=-10:1:20;                            %����ȷ�Χ
for SNR=-0:5:30                        
    sumtheta_1=0;                        %���
    for i=1:100                            
        xt=awgn(y1,SNR);                  
        [theta_1]=im_music(y1,SNR,length(DOA),M,N,fc,d,c); %music����
        theta_1=sort(theta_1); 
        sumtheta_1=sumtheta_1+(((0)-theta_1(1,1))^2+(10-theta_1(1,2))^2+(30-theta_1(1,3))^2+(50-theta_1(1,4))^2);
    end      
    RMS_1(k)=sqrt(sumtheta_1/4/100);   
    k=k+1;
end


%case2
ss1=[s1;s1;s2;conj(s2);s2;conj(s2)];
s_1=ss1(1:6,:);
DOA=[0,10,30,0,50,0];%�ź�����Ƕ�
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
snr=15;%�����
y2=awgn(y2,snr);%�������SNR������y����Ӧ�ĸ�˹������
R=y2*y2'/N;%Э�������

%case3
ss1=[s1;s1;s1;s2;conj(s2)];
s_1=ss1(1:5,:);
DOA=[0,10,30,50,0];%�ź�����Ƕ�
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
snr=15;%�����
y3=awgn(y3,snr);%�������SNR������y����Ӧ�ĸ�˹������
R=y3*y3'/N;%Э�������

%��Ӧ����
for k=1:p
    a1=ula(M,d,DOA(k),fc);
    A(:,k)=a1;
end

for k=1:p
    a1=ula(M,d,DOA(k),fc);
    A(:,k)=a1;
end

y1=[A*s;A'*s'];
snr=15;%�����
y1=awgn(y1,snr);%�������SNR������y����Ӧ�ĸ�˹������
R=y1*y1'/N;%Э�������

%��Բ�źŹ���
v=[A zeros(M,1);zeros(M,1) A'];
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