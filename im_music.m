function [DOA_guiji]= im_music(y,snr,p ,m,nn,fc,d,c)
y=awgn(y,snr);    %按规定信噪比加入噪声
Rxx=y*y'/nn;      
[s1,h1]=eig(Rxx);       %求出Rxx的特征向量与特征值
%%%估计噪声子空间%%%
if h1(1,1)>h1(2,2)      
    Vn=s1(:,p+1:2*m);     
else
    Vn=s1(:,1:2*m-p);            
end

%计算music谱
do=-90:90;
pu=zeros(1,length(do));
kg=1;
kt=2;
for k=-90:90
    a=zeros(2*m,2);
    for kk=1:m
        a(kk,1)=exp(-j*2*pi*fc*(kk-1)*d*sin(k/180*pi)/c);
        a(kk,2)=0;
    end
     for kk=m+1:2*m
       a(kk,1)=0;
       a(kk,2)=conj(exp(-j*2*pi*fc*(kk-m-1)*d*sin(k/180*pi)/c));
    end
    pu(1,kg)=1/det(a'*Vn*Vn'*a);%谱线
    kg=kg+1;
end

%侧向结果
for k=1:p
    [k1,k2]=max(pu);
    DOA_guiji(k)=(k2-1)-90;
    pu(k2)=0;
end