
function [yhat]=ilrma(x,fs,Ndft,K,L,Maxiter);
nmic=size(x,2);
M=floor(length(x(:,2))/250);
win=hamming(Ndft,'periodic');
OL=floor(Ndft*0.75);
xtf=[];
for i=1:nmic
[xtf(:,:,i),F]=STFTp(x(:,i),Ndft,win,OL,fs);
end
N=size(xtf,2);
Rx=[];
yhat=[];
delta=[];
Ka=[];
Wa=[];
Rphi=[];
epsi=1e-6;
phi=[];
z=[];
E=[];
T=[];
V=[];
g=[];
W=[];
D=[];
d=[];
U=[];
xtff=[];
Xp=[];
Wp=[];
ytf=[];
yhat=[];
P=[];
R=[];
Y=[];
lambda=[];
for f=1:length(F)
xtff{f}=squeeze(xtf(f,:,:)).';
Rx{f}=(1/N)*xtff{f}*xtff{f}';
[E{f},D{f}]=eig(Rx{f});
d{f} = real(diag(D{f}));
  [d{f} , ind] = sort(-d{f});
  E{f}=E{f}(:,ind(1:L));
  D{f}=diag(abs(d{f}(1:L)).^(-1/2));
  xtff{f}=xtff{f}-mean(xtff{f},2);
  Xp{f}=D{f}*E{f}'*xtff{f};
  Wp{f}=eye(L);
end

for f=1:length(F)
    ytf{f}=Wp{f}*Xp{f};
    Y(f,:,:)= abs(ytf{f}.').^2;
end 
for l=1:L
 P{l}=Y(:,:,l);
 T{l}=rand(length(F),K);
 V{l}=rand(K,N);
 R{l}=T{l}*V{l};
end

for iter=1:Maxiter
   
 for l=1:L
  T{l}=T{l}.*((((P{l}.*R{l}.^(-2))*V{l}.')./((R{l}).^(-1)*V{l}.')).^(1/2));
  R{l}=T{l}*V{l};
  V{l}=V{l}.*(((T{l}.')*((P{l}.*R{l}.^(-2)))./((T{l}.')*(R{l}).^(-1))).^(1/2)); 
  R{l}=T{l}*V{l};
  for f=1:length(F)
    U{l,f}=(Xp{f}*(Xp{f}'.*(((R{l}(f,:).').^(-1))*ones(1,L))))/N;
    truc=eye(L);
    Wa{f,l}=inv(Wp{f}*U{l,f})*truc(:,l);
     Wa{f,l}=Wa{f,l}*(Wa{f,l}'*U{l,f}*Wa{f,l})^(-1/2);
    
 end   
 end


  for f=1:length(F)
       Wp{f}=[];
          for l=1:L
         Wp{f}=[Wp{f} Wa{f,l}];
     end 
  end
  for f=1:length(F)
 ytf{f}=Wp{f}*Xp{f};
  Y(f,:,:)= abs(ytf{f}.').^2;
  end
  
  for l=1:L
    P{l}=Y(:,:,l);
  end
  
  for l=1:L
  lambda(l)=sqrt((1/(length(L)*length(F)))*sum(P{l},'all'));
  for f=1:length(F)
      Wp{f}(:,l)=Wp{f}(:,l)/lambda(l);
  end
  P{l}=P{l}/(lambda(l)^2);
  R{l}=R{l}/(lambda(l)^2);
  T{l}=T{l}/(lambda(l)^2);
  end 
 
end

for f=1:length(F)
    W{f}=Wp{f}*D{f}*E{f}';
    Ka{f}=pinv(W{f});
    W{f}=Ka{f}(1,1)*W{f};
    ytf{f}=W{f}*xtff{f};
    z(f,:,:)=ytf{f}.';
end
for i=1:L
[yhat(:,i)]=ISTFTp(z(:,:,i),size(x,1),win,OL);
end

