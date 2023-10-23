close all;
clear;
clc
%%
rng(130,'twister');
%%
addpath(genpath(fileparts(pwd)));
%%
M                   = 12;
Array               = 0:M-1;
KSamples            = 32;
SNR_dB              = 20; %[dB]
SNR_Lin             = 1/10^(SNR_dB / 20);
SIR_dB              = -10; %[dB]
SIR_Lin             = 1/10^(SIR_dB / 20);
ThetaVec            = linspace(0,180,1e3+1);
%%
ThetaDesired        = 80; %[deg]
ThetaInterf_1       = 50;
ThetaInterf_2       = 150;
%%

SigDesired          = randn(KSamples,1)+1j*randn(KSamples,1);
SigInterf_1Tmp      = randn(round(KSamples/2),1) + 1j*randn(round(KSamples/2),1);
SigInterf_2Tmp      = randn(round(KSamples/2),1) + 1j*randn(round(KSamples/2),1);
SigInterf_1         = SIR_Lin*[zeros(KSamples-length(SigInterf_1Tmp),1) ; SigInterf_1Tmp];
SigInterf_2         = SIR_Lin*[SigInterf_2Tmp ; zeros(KSamples-length(SigInterf_2Tmp),1)];
%%
SteeringVecDesired 	= exp(2j*pi*Array/2*cos(ThetaDesired/180*pi));
SteeringVecInterf_1 = exp(2j*pi*Array/2*cos(ThetaInterf_1/180*pi));
SteeringVecInterf_2 = exp(2j*pi*Array/2*cos(ThetaInterf_2/180*pi));

%%
SigDesiredAtArr            = transpose(SigDesired*SteeringVecDesired);
SigInterfAtArr_1           = transpose(SigInterf_1*SteeringVecInterf_1);
SigInterfAtArr_2           = transpose(SigInterf_2*SteeringVecInterf_2);

NoiseVec            = randn(M,KSamples);
KFactNoise         	= norm(NoiseVec(:,1));
KFactSig           	= norm(SigDesiredAtArr(:,1)); % For simplicity, here, the SNR is computed as follows. 
NoiseVec           	= SNR_Lin*NoiseVec / KFactNoise * KFactSig; %%

SigDesiredAtArr            = SigDesiredAtArr + SigInterfAtArr_1 + SigInterfAtArr_2 + NoiseVec;
%%
SigSeg_1            = SigDesiredAtArr(:,1:end/2);
SigSeg_2            = SigDesiredAtArr(:,end/2:end);
CorrSeg_1           = 1/KSamples * SigSeg_1*SigSeg_1';
CorrSeg_2           = 1/KSamples * SigSeg_2*SigSeg_2';
CorrTensor(:,:,1)   = CorrSeg_1;
CorrTensor(:,:,2)   = CorrSeg_2;
CorrEuc             = mean(CorrTensor,3);
CorrRiem            = RiemannianMean(CorrTensor); 
%%
for i = 1:length(ThetaVec)
    SteeringVecSpectrum        = exp(2j*pi*Array'/2*cos(ThetaVec(i)/180*pi));
    MLSpectrumEuc(i)           = SteeringVecSpectrum'*CorrEuc*SteeringVecSpectrum;
    MLSpectrumRiem(i)          = SteeringVecSpectrum'*CorrRiem*SteeringVecSpectrum;   
end

%%
MLSpectrumEucNorm = NormSpectrumFunc(MLSpectrumEuc);
MLSpectrumEucRiem = NormSpectrumFunc(MLSpectrumRiem);

rLimScale = [-20 0 ];

figure;
polarplot(ThetaVec/180*pi,MLSpectrumEucRiem,'LineWidth',3); hold on%
s = polarplot(ThetaVec/180*pi,MLSpectrumEucNorm,':','LineWidth',3);
s.Color(4)=0.9;
p = polarplot(ThetaDesired/180*pi*ones(1,2),rLimScale,'LineWidth',2);
p.Color = 'black';
p = polarplot(ThetaInterf_1/180*pi*ones(1,2),rLimScale,'--','LineWidth',2);
p.Color = 'black';
p = polarplot(ThetaInterf_2/180*pi*ones(1,2),rLimScale,'--','LineWidth',2);
p.Color = 'black';
rlim(rLimScale);
%
legend('\Gamma_R (Riem)','\Gamma_E (Euc)','Desired','Interference');
thetalim([0 180])


%%
function SpectrumFull = NormSpectrumFunc(Spectrum)
SpectrumFull = 10*log10(abs(Spectrum)/max(abs(Spectrum)));
end


