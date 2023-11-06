[x,fs]=audioread('PFE test.wav');% Reading the mixture
x=10*x;
sigma=0.000001;
x=x+sigma^2*randn(size(x,1),size(x,2));
mic=size(x,2); % Finding the number of microphones
T=size(x,1); % Finding the number of samples
eta0=3; % Initializing the leraning rate
L=2; % Setting up the number of estimated sources
Ndft=512; % Initializing the number of fast fourier transform's points
OL=floor(0.5*Ndft); % Initializing the overlap size
win=2*hann(Ndft,'periodic')/Ndft; % Initializing the window
SHIFT = Ndft - OL; 
cmp=1; %Initializing the time counter
beta=0.5; %Initializing the forgetting factor
Y = zeros(T+Ndft+Ndft, L); %Initializing the variable which will contain the output separated sources
Wa = zeros(T+Ndft+Ndft, L); %Initializing the window of inverse STFT
moy=zeros(mic,Ndft/2+1); % Initializing the matrix which wil contain the mean for each vector x through time
norma=zeros(L,1);


%Initializing the matrix for separation
for i = 1:Ndft/2 +1
    Rx{i}=zeros(mic,mic);
    Rphi{i}=zeros(L,L);
    Wp{i}=eye(L);
    Ksi(i)=0.0000001;
end

N=(fix((size(x,1)+Ndft)/SHIFT)-1); % Giving the number of time samples

% Beginning the adaptive algorithm through time
while cmp<N+1


  % computing the STFT of each vector x_m
  for i=1:mic

        X=x(:,i);
        Xautre= [zeros(Ndft, 1); X; zeros(Ndft, 1)]; %zeros padding
        sp = SHIFT*cmp + 1; %computing the position of shifting in order to scroll the window
        Xn = win.*Xautre(sp:sp+Ndft-1); %multiplying by the window
        xtf{i}= fft(Xn); %computing the FFT of the resulting signal
        xtf{i} = xtf{i}(1:fix(Ndft/2)+1); %Taking only the positive frequencies 
        F = [0:fix(Ndft/2)]'.*fs/Ndft; %computing the number of frequencies
  end

    xtff=cell2mat(xtf).'; % Reshaping the data into a matrix mic-by-F
 
    moy=((cmp-1)/cmp)*moy+xtff/cmp; % computing the mean 

    % Whitenning the data
    if cmp==1
     xtfff=xtff;
    else
    xtfff=xtff-moy;
    end
  for f=1:length(F)

    Rx{f}=((cmp-1)/cmp)*Rx{f}+xtfff(:,f)*xtfff(:,f)'/cmp;
    [E{f},D{f}]=eig(Rx{f});
    d = real(diag(D{f}));
    [d , ind] = sort(-d);
    E{f}=E{f}(:,ind(1:L));
    D{f}=diag(abs(d(1:L)).^(-1/2));
  %D{f}=D{f}(:,ind(1:L));
                            
    Xp(:,f)=D{f}*E{f}'*xtfff(:,f);
    
 
 
    ytf(:,f)=Wp{f}*Xp(:,f);
  end  

% Computing the estimated sources norm 
for l=1:L
    norma(l)=norm(ytf(l,:));
end
% Beginning the minimizing of the natural gradient
for f=1:length(F)
     phi{f}=rdivide(ytf(:,f),norma+(1e-6)*ones(L,1)); %computing the score function
     Rphi{f}=phi{f}*ytf(:,f)'; %computing the instantaneous covariance matrix
     Delta{f}=((diag(diag(Rphi{f}))-Rphi{f})*Wp{f}); %computing the natural gradient
     Ksi(f)=beta*Ksi(f)+(1-beta)*norm(Xp(:,f))^2/L; %normalizing the learning rate
     Wp{f}=Wp{f}+eta0*sqrt(1/(Ksi(f)+1e-6))*Delta{f};
 
end

 % Rescaling
 for f=1:length(F)

    W{f}=Wp{f}*D{f}*E{f}';
    Ka{f}=pinv(W{f});
    W{f}=diag(diag(Ka{f}))*W{f};

    ytf(:,f)=W{f}*xtff(:,f);

end
% Computing the STFT inverse
for i=1:L

Yautre=[ytf(i,:).' ;conj(ytf(i,end-1:-1:2).')]; % Rebuild the negative frequencies
tmp = real(ifft(Yautre)); % Computing the FFT inverse
Wa(sp:sp+Ndft-1,i) = Wa(sp:sp+Ndft-1,i) + win.^2; % Computing the window for inverse STFT
Y(sp:sp+Ndft-1,i) = Y(sp:sp+Ndft-1,i) + win.*tmp(1:Ndft); % Rebuild our separated signal by taking the mean in case of overlap

end




cmp=cmp+1;      
end

for i=1:L
    Wa(SHIFT*(N+1)+1:SHIFT*(N+1)+Ndft,i) = Wa(SHIFT*(N+1)+1:SHIFT*(N+1)+Ndft,i) + win.^2;
    Yhat(:,i) = Y(Ndft+1:Ndft+T,i)./Wa(Ndft+1:Ndft+T,i);
end 
for i=1:L
    filename = sprintf('source %d.wav', i);
    audiowrite(filename, Yhat(:, i), 16000);
end

