[x,fs]=audioread('PFE test.wav');
x=10*x(:,[4 10]);
sigma=0.0001;
eta1=0.000001;
eta2=0.000001;
mic=size(x,2);
alpha=0.9;
x=x+sigma^2*randn(size(x,1),mic);
T=size(x,1);
L=2;
Ndft=4096 ;
OL=floor(0.75*Ndft);
win=hann(Ndft,'periodic');
xtf=[];
SHIFT = Ndft - OL;
cmp=1;
beta=0.9;
norma=zeros(L,1);
Xp=[];
epsi=0.0000001;
ytff=[];
ytf=[];
yhat=[];
Rphi=[];
g=[];
Ksi=[];
Y = zeros(T+Ndft+Ndft, L);
Wa = zeros(T+Ndft+Ndft, L);
moy=zeros(mic,Ndft/2 +1);

for i = 1:Ndft/2 +1
    Rx{i} = 0.001*eye(mic,mic);
    Rphit{i}=epsi*eye(L);
    Wp{i}=eye(L);

end
for l=1:L
for f= 1:Ndft/2 +1
    U{l,f}=(1e-9)*eye(L);
end
end

E=eye(L);
N=(fix((size(x,1)+Ndft)/SHIFT)-1);

sp=1;
while cmp<N+1
   
  for i=1:mic
        X=x(:,i);
        Xautre=[zeros(Ndft,1); X; zeros(Ndft,1)];
        sp = SHIFT*cmp + 1;
        Xn = win.*Xautre(sp:sp+Ndft-1);
        xtf{i}= fft(Xn);
        xtf{i} = xtf{i}(1:fix(Ndft/2)+1);
        F = [0:fix(Ndft/2)]'.*fs/Ndft;
  end
    xtfff=cell2mat(xtf).';
    
 
 
 for iter=1:2
   r=ones(L,1)*epsi;
    for l=1:L
     for f=1:length(F)
    r(l)=r(l)+ abs(Wp{f}(l,:)*xtfff(:,f))^2;
     end
     r(l)=sqrt(r(l));
     for f=1:length(F)
    U{l,f}=alpha*U{l,f}+((1-alpha)/(2*r(l)+epsi))*xtfff(:,f)*xtfff(:,f)';
     end
    end
    
     for l=1:L
     for m=1:mic
         if l==m
           for f=1:length(F)
        V(m,l,f)=1-(Wp{f}(l,:)*U{l,f}*Wp{f}(l,:)'+epsi)^(-0.5);
           end
         else
           for f=1:length(F)
             V(m,l,f)=(Wp{f}(m,:)*U{m,f}*Wp{f}(l,:)')/(Wp{f}(l,:)*U{m,f}*Wp{f}(l,:)'+epsi);
           end
         end
     end
 for f=1:length(F)
    Wp{f}=Wp{f}-V(:,l,f)*Wp{f}(l,:);
 end
     end
 end
 
 for f=1:length(F)
    W{f}=Wp{f};
    Ka{f}=inv(W{f});
    W{f}=diag(diag(Ka{f}))*W{f};
    ytff(:,f)=W{f}*xtfff(:,f);
end
 

for i=1:L
Yautre=[ytff(i,:).' ;conj(ytff(i,end-1:-1:2).')];
tmp = real(ifft(Yautre));
Wa(sp:sp+Ndft-1,i) = Wa(sp:sp+Ndft-1,i) + win.^2;
Y(sp:sp+Ndft-1,i) = Y(sp:sp+Ndft-1,i) + win.*tmp(1:Ndft);

end
 


cmp=cmp+1;      
end

for i=1:L
    Wa(SHIFT*(N+1)+1:SHIFT*(N+1)+Ndft,i) = Wa(SHIFT*(N+1)+1:SHIFT*(N+1)+Ndft,i) + win.^2;
    Yjdid(:,i) = Y(Ndft+1:Ndft+T,i)./Wa(Ndft+1:Ndft+T,i);
end
