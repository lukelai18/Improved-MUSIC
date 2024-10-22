%����QPSK�ź�
function [data_c]=qpsk(num,f)
data=rand(1,num);
%figure(1)
%plot(data)
%title('����ʱ����');
Rb=20000;
Ts=1/f;
Ns=1000;
sample=1*Ns;
N=sample*length(data)/2;
data1=2*data-1;

data_1=zeros(1,N);
for i=1:num/2
    data_1(sample*(i-1)+1:sample*i)=data1(2*i-1);
end
data_2=zeros(1,N);
for i=1:num/2
    data_2(sample*(i-1)+1:sample*i)=data1(2*i);
end
a=zeros(1,N);
b=zeros(1,N);
for i=1:N
    a(i)=cos(2*pi*f*(i-1)*Ts/Ns);
    b(i)=sin(2*pi*f*(i-1)*Ts/Ns);
end
data_a=data_1.*a;
data_b=data_2.*b;
data_c=data_a+j*data_b;
%figure(2)
%subplot(311)
%plot(data_a)
%title('QPSK�ѵ�ʵ��ʱ���ź�');
%subplot(312)
%plot(data_b)
%title('QPSK�ѵ��鲿ʱ���ź�');
%subplot(313)
%plot(data_c)
%title('QPSK�ѵ��ź�ʱ����');