%���þ������о���
function av=ula(m,d,theta,fc)
%mΪ��Ԫ������dΪ��Ԫ��࣬thetaΪ�ź�����Ƕȣ�fcΪ��Ƶ
%theta=30;m=3;fc=500;  %������Ԫ���Ϊlamda/2������Ƕ�Ϊ30�㣬��ƵΪ500Hz,��Ԫ����Ϊ3
%���ô������ľ���
%b=zeros(m,2);
%for n=1:m
%b(n,1)=(n-1)*100;
%b(n,2)=0;
%end
lamda=3e8/fc;%����
%av=zeros(m,1);
a=exp(-1i*2*pi*d*sind(theta)/lamda);
for n=1:m
    %av(n)=exp(j*2*pi*(b(n,1)*sind(theta)+b(n,2)*cosd(theta))/lamda);
    av(n)=a.^(n-1);
end
av=av';
