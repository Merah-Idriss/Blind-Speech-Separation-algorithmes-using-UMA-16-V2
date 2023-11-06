% Initializing metrics vectors
SIR=[];
SAR=[];
SDR=[];
sirt=[];
sdrt=[];
sart=[];
permt=[];
PERM=[];
%Normalizing the sources
[s1,fs]=audioread('Dris.wav');
s1=s1/norm(s1);
[s2,fs]=audioread('Zakariya.wav');
s2=s2/norm(s2);
se=Yjdid;
Te=min([size(s1,1) size(s2,1) size(se,1)]);
s=[s1(1:Te,:),s2(1:Te,:)].';
se=se(1:Te,:).';
se(1,:)=se(1,:)/norm(se(1,:));
se(2,:)=se(2,:)/norm(se(2,:));
%Computing the metrics over time
sp=1;
pas=27;
inc=floor(Te/pas);
for i=1:pas
    [sirt,sart,sdrt,permt]=bss_eval_sources(se(:,sp:(sp+inc-1)),s(:,sp:(sp+inc-1)));
    SIR=[SIR sirt];
    SDR=[SDR sdrt];
    SAR=[SAR sart];
    PERM=[PERM permt];
    sp=sp+inc;
end
% Computing the mean after convergence 
SARmoy=mean(SAR(:,5:end),2);
SDRmoy=mean(SDR(:,5:end),2);
SIRmoy=mean(SIR(:,5:end),2);
