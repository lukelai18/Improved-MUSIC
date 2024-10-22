%…Ë÷√BPSK–≈∫≈
function [bpsk]=bkskdigital(s,f)
t=0:0.4*pi/99:2*pi+1.6*pi/99;
l=length(t);
cp=[];
mod=[];mod1=[];bit=[];

for n=1:length(s)
    if s(n)==0
        cp1=-ones(1,500);
        bit1=zeros(1,500);
    else s(n)==1
         cp1=ones(1,500);
        bit1=ones(1,500);
    end
    c=sin(f*t);
    cp=[cp cp1];
    mod=[mod c];
    bit=[bit bit1];
end
bpsk=cp.*mod;
 
