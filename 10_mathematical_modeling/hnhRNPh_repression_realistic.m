clear
% close all

% This code implements a three-tiered cascade, in which each level is
% described by an algebraic equation (see Methods and Supplemental Material). 
% It calculates dose-response curves and Hill coefficients for an example 
% parameter set (i=1) and for randomly sampled parameter values

%-------------------------

% 1) Definition of parameter values

% Cascade parameters
K10=1;
K20=0.25; % 0.05
nH10=2.2;
nH20=2;

alpha0=0.01; % kon repression factor by hnRNPh
ksplAE0=3; % splicing rate of the AE
kCE0=1; %Splicing rate of the constitutive exon
tau0=7; % Time delay from co-transcriptional splicing

% % Cascade parameters
% K10=1;
% K20=0.25; % 0.05
% nH10=2.2;
% nH20=2;
% 
% alpha0=0.01; % kon repression factor by hnRNPh
% ksplAE0=3; % splicing rate of the AE
% kCE0=1; %Splicing rate of the constitutive exon
% tau0=7; % Time delay from co-transcriptional splicing

% Parameters for Sampling
mu=0;
sig=0.73;
signH=0.2;

%-------------------------

% 2) Input titration and parameter sampling

for i=1:2 % Example parameter set (i=1) or random parameter sampling (i=2)

if i==1 % Example simulation
   
    freeRBP=logspace(-10,10,5000); % titration of HnRNPh levels
    K1=K10; % reference parameter values
    K2=K20;
    alpha=alpha0;
    ksplAE=ksplAE0;
    kCE=kCE0;
    tau=tau0;

    nH1=nH10;
    nH2=nH20;


elseif i==2 % Random parameter sampling
    
    figure
    freeRBP=logspace(-10,10,5000);  % titration of HnRNPh levels
    no=1000; % Sample size
    
     % Sampling the half-maximal input at each level 
    K1=K10.*lognrnd(mu,sig,no,1);
    K2=K20.*lognrnd(mu,sig,no,1);

    alpha=alpha0.*lognrnd(mu,sig,no,1);
    ksplAE=ksplAE0.*lognrnd(mu,sig,no,1);
    kCE=kCE0.*lognrnd(mu,sig,no,1);
    tau=tau0.*lognrnd(mu,sig,no,1);

    % Sampling the Hill coefficient at each level 
    nH1=0.4*randn(no,1)+nH10;
    nH2=0.4*randn(no,1)+nH20;

end

%---------------------

% 3) Equations describing the three-tiered cascade
% input: hnRNPh concentration ('freeRBP'); output: PSI

boundRBP=freeRBP.^nH1./(freeRBP.^nH1+K1.^nH1); % Level 1: cooperative hnRNPh binding

exondef=(1+(alpha-1).*boundRBP).^nH2./((1+(alpha-1).*boundRBP).^nH2+K2); % Level 2: Exon definition controlled by hnRNPh

kAE=ksplAE.*exondef;
PSI=(1-(kCE./(kCE+kAE)).*(exp(-kAE.*tau))); % Level 3: PSI controlled by exon definition status

%----------------------

if i==1 % Specific parameter set

% 4) Calculate Hill coefficients at each level and globally for the cascade
% using the formula nH = log81/log(EC90/EC10)
% EC10 and EC90: inputs required for 10% and 90% activation

% Local Hill coefficient: Level 1
EC10=freeRBP(min(find(boundRBP>0.1)));
EC90=freeRBP(min(find(boundRBP>0.9)));
nH=log(81)./log(EC90./EC10);

% Local Hill coefficient: Level 2
EC10a=boundRBP(min(find((exondef-min(exondef))./(max(exondef)-min(exondef))<0.1)));
EC90a=boundRBP(min(find((exondef-min(exondef))./(max(exondef)-min(exondef))<0.9)));
nHa=log(81)./log(EC90a./EC10a);

% Local Hill coefficient: Level 3
EC10b=exondef(min(find((PSI-min(PSI))./(max(PSI)-min(PSI))<0.1)));
EC90b=exondef(min(find((PSI-min(PSI))./(max(PSI)-min(PSI))<0.9)));
nHb=log(81)./log(EC90b./EC10b);

% Global Hill coefficient: change in RBPfree required to switch PSI
EC10c=freeRBP(min(find((PSI-min(PSI))./(max(PSI)-min(PSI))<0.1)));
EC90c=freeRBP(min(find((PSI-min(PSI))./(max(PSI)-min(PSI))<0.9)));
nHc=log(81)./log(EC90c./EC10c);

%-------------------------------

% 5) Plotting the dose-response curves

% Level 1
subplot(2,2,1)
semilogx(freeRBP,boundRBP,'LineWidth',3)
title(['step 1; nH=' num2str(nH)])
xlabel('[hnRNPh]')
ylabel('bound')
xlim([0.1*EC10 10*EC90])
ylim([0 1])
hold on

% Level 2
subplot(2,2,2)
plot(boundRBP,exondef,'LineWidth',3)
title(['step 2; nH=' num2str(nHa)])
xlabel('bound')
ylabel('%exon defined')
% xlim([0.1*EC90a 10*EC10a])
ylim([0 1])
hold on

% Level 3
subplot(2,2,3)
plot(exondef,PSI,'LineWidth',3)
title(['step 3; nH=' num2str(nHb)])
xlabel('%exon defined')
ylabel('PSI')
hold on
% xlim([0.1*EC10b 10*EC90b])

% Overall cascade
subplot(2,2,4)
semilogx(freeRBP,PSI,'LineWidth',3)
title(['overall; nH=' num2str(nHc)])
xlabel('[hnRNPh]')
ylabel('PSI')
hold on
xlim([0.1*EC90c 10*EC10c])

%----------------------

elseif i==2 % Parameter sampling

% 4a) Calculate Hill coefficients at each level and globally for the cascade
% 5a) Plot overall dose-response curve of the cascade

% Plot non-normalized dose-response curves
subplot(1,3,1),semilogx(freeRBP,PSI,'LineWidth',2)
xlabel('[hnRNPh]')
ylabel('PSI')
xlim([0.01*K10 1000*K10])
ylim([0 1])
hold on
title('Cascade dose-response')

% Plot min-max normalized dose-response curves
PSI_norm=(PSI'-repmat(min(PSI'),size(PSI',1),1))./(repmat(max(PSI'),size(PSI',1),1)-repmat(min(PSI'),size(PSI',1),1));
subplot(1,3,2),semilogx(freeRBP,PSI_norm')
xlim([0.01*K10 1000*K10])
hold on
[x10,y10]=min(abs(PSI_norm-0.1));
[x90,y90]=min(abs(PSI_norm-0.9));
xlabel('[hnRNPh]')
ylabel('min-max normalized PSI')
title('Normalized cascade dose-response')

% Calculate Hill coefficients and plot distribution
EC10=freeRBP(y10);
EC90=freeRBP(y90);
subplot(1,3,3),hist(log(81)./log(EC90./EC10),no/100)
xlabel('nH')
ylabel('%')
title(['nH distribution: median=' num2str(median(log(81)./log(EC90./EC10)))])

end

end

