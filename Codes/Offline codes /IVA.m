function y=IVA(x,fs,eta,Ndft,L,activation,Maxiter);
% IVA offline blind source separation 
% 
% Sources are estimated from their mixtures using the natural gradient descent
%
% [y]=IVA(x,fs,eta,Ndft,L,activation)
%
% Inputs:
% x: nmic x nsampl matrix containing mixtures
% fs: sampling frequency Hz
% eta: learning rate for the gradient descent <1
% Ndft: number of frequency bins for the STFT pow(2)
% activation : the score function of the sources 

% Outputs:
% y: nsrc x nsampl matrix containing estimated signals


nmic=size(x,2);
OL=floor(0.75*Ndft);
win=hann(Ndft,'periodic');
xtf=[];

% Compute the Short-Time Fourier Transform for each channel using STFTp 
% it returns xtf (time frequency domaine tensor F x N x nmic
for i=1:nmic
[xtf(:,:,i),F]=STFTp(x(:,i),Ndft,win,OL,fs);
end
% N is the number of time frames
N=size(xtf,2);
% Initialize the matrices
Rx=[];
yhat=[];
delta=[];
Ka=[];
Rphi=[];
epsi=1e-6;
p=[];
phi=[];
z=[];
E=[];
g=[];
W=[];
D=[];
d=[];
xtff=[];
Xp=[];
Wp=[];
ytf=[];
yhat=[];

% Perform whitenning for each frequency bin f
for f=1:length(F)
xtff{f}=squeeze(xtf(f,:,:)).';
Rx{f}=(1/(N))*xtff{f}*xtff{f}';
[E{f},D{f}]=eig(Rx{f});
d{f} = real(diag(D{f}));
  [d{f} , ind] = sort(-d{f});
  E{f}=E{f}(:,ind(1:L));
  D{f}=diag(abs(d{f}(1:L)).^(-1/2));
  xtff{f}=xtff{f}-mean(xtff{f},2);
  Xp{f}=D{f}*E{f}'*xtff{f};
  Wp{f}=eye(L);
end

%choose activation function
for iter=1:Maxiter
if activation=='regular'
 %compute the output ytf
    for f=1:length(F)
        ytf{f}=Wp{f}*Xp{f};
        g(:,:,f)=ytf{f};
    end
    %compute the norms matrix norma (nsrc x Nframes) 
    norma=zeros([L,N]);
    for n=1:N
     for l=1:L
       for f=1:length(F)
        norma(l,n)=norma(l,n)+abs(g(l,n,f))^2;   
       end
       norma(l,n)=sqrt(norma(l,n));  
        norma(l,n)=sqrt(norma(l,n)+epsi);  
     end         
    end
phi=[];

%use the gradient descent to update Wp for each f
 for f=1:length(F)
    %calculate the score matrix for (nsrc x Nframes)
     phi{f}=rdivide(ytf{f},norma); 
     %estimate the correlation across all time frames
     Rphi{f}=(phi{f}*ytf{f}')/N;
     %update Wp 
     Wp{f}=Wp{f}+eta*((eye(L)-Rphi{f})*Wp{f});
 end
 
elseif activation=='sigmoid'
     for f=1:length(F)
        ytf{f}=Wp{f}*Xp{f};
    end
  
phi=[];
 for f=1:length(F)
     phi{f}=2*Sigmoid(ytf{f})-ones(L,N); 
     Rphi{f}=(phi{f}*ytf{f}')/N;
     Wp{f}=Wp{f}+eta*((eye(L)-Rphi{f})*Wp{f});
 end
end
  

end
% remove whitenning and fix scaling ambiguity
for f=1:length(F)
    W{f}=Wp{f}*D{f}*E{f}';
    Ka{f}=pinv(W{f});
    W{f}=diag(diag(Ka{f}))*W{f};
    ytf{f}=W{f}*xtff{f};
    z(f,:,:)=ytf{f}.';
end
%compute ISTFT to obtain the estimated sources
for i=1:L
[y(:,i)]=ISTFTp(z(:,:,i),size(x,1),win,OL);
end


