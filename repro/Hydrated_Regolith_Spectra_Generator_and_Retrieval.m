%%
clear all
close all
clc
%%
%Replace with location of spectra and script files
%cd 'C:\Users...'
%% Define wavelength range and spacing for simulations
WLS=linspace(1,4,601).';
%% Input the total number of simulations to run. first 1/2 will be mare, second 1/2 will be highlands
NumSpectra=100;
%% Set abundance limits for spectral endmembers in mixture
LowAbunReg=0.0;
HighAbunReg=0.25;

LowAbunPyr=0.0;
HighAbunPyr=0.25;

LowAbunMORB=0.0;
HighAbunMORB=0.3;

rng(0,'twister');
%% Load in optical constants or reflectance data and calculate SSA using Hapke model

%Density and grain size of MORB glass
densMORB=2.8;
dMORB=69E-6;

%Low wavelength portion of MORB spectrum (<1.5 microns)
Morb_D38A_LowLam=readtable('../data/Morb_D38A_Low_wavelength.txt');
Morb_D38A_LowLam=table2array(Morb_D38A_LowLam);

WaterPPM=[1522,762,176,22]; %Total water measured from step-heating experiments

%Load in MORB step-wise heating spectra
MORB_D38A=dir('../data/*.csv');
%Use the 650C spectrum for the lower wavelength portion of the spectrum
MORB_D38A2=fullfile(MORB_D38A(1).folder, MORB_D38A(1).name);
MORB_D38A2=readtable(MORB_D38A2);
MORB_D38A2=table2array(MORB_D38A2);
MORB_D38A2(:,1)=1E4./MORB_D38A2(:,1); %convert from cm^-1 to microns
MORB_D38A2(:,2)=MORB_D38A2(:,2)/100; %convert from percent to decimal
MORB_D38A2=MORB_D38A2(6481:11648,:); %only use the spectral region of interest
MORB_D38A2=flipud(MORB_D38A2);
MORB_D38A_MidLam=MORB_D38A2;

%Normalize all MORB spectra to the reflectance at 2.6 microns at the
%lowest water amount and highest temperature (22 ppm, 800 C)
NormalizationValue=0.2363;
NormFactor=MORB_D38A_MidLam(end,2)/NormalizationValue;
MORB_D38A_MidLam(:,2)=MORB_D38A_MidLam(:,2)./NormFactor;
clear MORB_D38A2
%Perform Normalization of MORB spectra at 650, 700, 750, 800 C
for i=1:4
    MORB_D38A2=fullfile(MORB_D38A(i).folder, MORB_D38A(i).name); %start from 650 C spectrum
    MORB_D38A2=readtable(MORB_D38A2);
    MORB_D38A2=table2array(MORB_D38A2);
    MORB_D38A2(:,1)=1E4./MORB_D38A2(:,1); %convert from cm^-1 to microns
    MORB_D38A2(:,2)=MORB_D38A2(:,2)/100; %convert from percent to decimal
    MORB_D38A2=flipud(MORB_D38A2);
    %Normalize
    NormFactor=MORB_D38A2(12752,2)/NormalizationValue; %normalize at 2.6 micron
    MORB_D38A2(:,2)=MORB_D38A2(:,2)./NormFactor;
    %Stitch with low wavelength spectrum
    Morb_D38A_LowLam2=Morb_D38A_LowLam;
    Morb_D38A_LowLam2(:,2)=Morb_D38A_LowLam(:,2)*(MORB_D38A2(6901,2)/0.15);
    Morb_D38A_LowLam2=Morb_D38A_LowLam2(1:15,:);
    MORB_D38A2=cat(1,Morb_D38A_LowLam2,MORB_D38A2(6902:end,:));
    % Interpolate data to 5 nm spacing between 1 and 4 microns
    InterpMorb=interp1(MORB_D38A2(:,1),MORB_D38A2(:,2),WLS);
    if i>1 % Use 650 C spectrum below 2.6 microns to isolate 3 micron feature changes
        InterpMorb(1:321)=RMORB(1:321,1);
    else
    end
    RMORB(:,i)=InterpMorb;
    %Use Hapke model to convert from laboratory reflectance to single
    %scattering albedo using reflectance and scattering asymmetry
    % factor (p) of 0.81 (see manuscript for details on Hapke model)
    SSAMORB(:,i)=Hapke_Inverse_Function_Passive(RMORB(:,i),0.81,WLS);
end
%% Interpolation of MORB Spectra
%Set range and resolution for interpolation
%(standard is 1 ppm spacing from 0 to 1666 ppm) With 30% maximum MORB
%abundance this gives a maximum total water abundance of 500 ppm.
WatInterp=linspace(0,1666,1667);

for m=1:numel(WLS)
    %For each wavelength, fit a line to the ESPAT values for each heating step spectrum
    ESPAT(m,:)=polyfit(WaterPPM(1:4),(1-SSAMORB(m,1:4))./SSAMORB(m,1:4),1); %linear function
    %For each wavelength, create new synthetic spectra from the linear ESPAT relationship
    TotalWatInterpESPAT(:,m)=1./(ESPAT(m,1).*WatInterp+ESPAT(m,2)+1);
end

% Check match of interpolated spectra with real MORB spectra
figure(1);
RealMorbIndices=[23,177,763,1523];
for k=1:4
    scatter(WLS,SSAMORB(:,k),20,'filled');
    hold on
end
for k=1:4
    plot(WLS,TotalWatInterpESPAT(RealMorbIndices(k),:),'Color','black','linestyle','-','linewidth',2);
end
xlabel('Wavelength (\mum)','FontSize',16,'FontWeight','bold');
ylabel('SSA','FontSize',16,'FontWeight','bold');
ax=gca;
ax.XLim=[2.5,3.5];
ax.LineWidth = 2;
ax.TickDir='in';
ax.FontName='arial';
ax.FontSize=14;
ax.FontWeight='bold';
box on
grid on
%% Load in other spectral endmembers (regolith and pyroxene)
densReg=1.8;
dReg=32E-6;
%Mare Mature Endmember
MareMature=readtable('../data/Mare_70181_Spectra.txt');
MareMature=table2array(MareMature);
MareMatureRaw=MareMature;
%Remove Water and organics Signature by fitting a linear continuum between
%2.65 and 3.8 microns
[V ind32]=min(abs(MareMature(:,1)-2.65));
[V ind36]=min(abs(MareMature(:,1)-3.8));
SandIQuad=polyfit([MareMature(ind32,1),MareMature(ind36,1)],[MareMature(ind32,2),MareMature(ind36,2)],1);
NewRs=MareMature(ind32:ind36,1)*SandIQuad(1)+SandIQuad(2);
MareMature(ind32:ind36,2)=NewRs;

PassiveMareMature=interp1(MareMature(:,1),MareMature(:,2),WLS);
SSAMareMature=Hapke_Inverse_Function_Passive(PassiveMareMature,0.81,WLS);

% Mare Immature Endmember
MareImmature=readtable('../data/Mare_71061_Spectra.txt');
MareImmature=table2array(MareImmature);
MareImmatureRaw=MareImmature;
%Remove Water and organics Signature
[V ind32]=min(abs(MareImmature(:,1)-2.65));
[V ind36]=min(abs(MareImmature(:,1)-3.8));
SandIQuad=polyfit([MareImmature(ind32,1),MareImmature(ind36,1)],[MareImmature(ind32,2),MareImmature(ind36,2)],1);
NewRs=MareImmature(ind32:ind36,1)*SandIQuad(1)+SandIQuad(2);
MareImmature(ind32:ind36,2)=NewRs;

PassiveMareImmature=interp1(MareImmature(:,1),MareImmature(:,2),WLS);
SSAMareImmature=Hapke_Inverse_Function_Passive(PassiveMareImmature,0.81,WLS);

% Highlands Mature Endmember
HighlandsMature=readtable('../data/Highlands_62231_Spectra.txt');
HighlandsMature=table2array(HighlandsMature);
HighlandsMatureRaw=HighlandsMature;
%Remove Water and organics Signature
[V ind32]=min(abs(HighlandsMature(:,1)-2.65));
[V ind36]=min(abs(HighlandsMature(:,1)-3.8));
SandIQuad=polyfit([HighlandsMature(ind32,1),HighlandsMature(ind36,1)],[HighlandsMature(ind32,2),HighlandsMature(ind36,2)],1);
NewRs=HighlandsMature(ind32:ind36,1)*SandIQuad(1)+SandIQuad(2);
HighlandsMature(ind32:ind36,2)=NewRs;

PassiveHighlandsMature=interp1(HighlandsMature(:,1),HighlandsMature(:,2),WLS);
SSAHighlandsMature=Hapke_Inverse_Function_Passive(PassiveHighlandsMature,0.81,WLS);

% Highlands Immature Endmember
HighlandsImmature=readtable('../data/Highlands_61221_Spectra.txt');
HighlandsImmature=table2array(HighlandsImmature);
HighlandsImmatureRaw=HighlandsImmature;
%Remove Water and organics Signature
[V ind32]=min(abs(HighlandsImmature(:,1)-2.65));
[V ind36]=min(abs(HighlandsImmature(:,1)-3.8));
SandIQuad=polyfit([HighlandsImmature(ind32,1),HighlandsImmature(ind36,1)],[HighlandsImmature(ind32,2),HighlandsImmature(ind36,2)],1);
NewRs=HighlandsImmature(ind32:ind36,1)*SandIQuad(1)+SandIQuad(2);
HighlandsImmature(ind32:ind36,2)=NewRs;

PassiveHighlandsImmature=interp1(HighlandsImmature(:,1),HighlandsImmature(:,2),WLS);
SSAHighlandsImmature=Hapke_Inverse_Function_Passive(PassiveHighlandsImmature,0.81,WLS);

%Apollo 15 Pyroxene Endmember
densPyrox=3.2;
dPyrox=70E-6;
Apollo15_Pyroxene=readtable('../data/Apollo15Sample15555ReddishBrownPyroxeneB.txt');
Apollo15_Pyroxene=table2array(Apollo15_Pyroxene);
Apollo15_PyroxeneRaw=Apollo15_Pyroxene;
%Remove Water and organics Signature
[V, ind32]=min(abs(Apollo15_Pyroxene(:,1)-3.2)); %only need to remove organics here
[V ind36]=min(abs(Apollo15_Pyroxene(:,1)-3.6));
SandIQuad=polyfit([Apollo15_Pyroxene(ind32,1),Apollo15_Pyroxene(ind36,1)],[Apollo15_Pyroxene(ind32,2),Apollo15_Pyroxene(ind36,2)],1);
NewRs=Apollo15_Pyroxene(ind32:ind36,1)*SandIQuad(1)+SandIQuad(2);
Apollo15_Pyroxene(ind32:ind36,2)=NewRs;

PassivePyrox=interp1(Apollo15_Pyroxene(:,1),Apollo15_Pyroxene(:,2),WLS);
RPyrox=PassivePyrox;
SSAPyrox=Hapke_Inverse_Function_Passive(RPyrox,0.81,WLS);
%%
%Place densities, grain sizes, and SSAs into arrays
rho=1000*[NaN,densMORB,densMORB,densMORB,densMORB,densReg,densReg,densReg,densReg,densPyrox];
d=[NaN,dMORB,dMORB,dMORB,dMORB,dReg,dReg,dReg,dReg,dPyrox];
ConstituentSSAs=[SSAMORB,SSAHighlandsMature,SSAHighlandsImmature,SSAMareMature,SSAMareImmature,SSAPyrox];
%% Create spectral mixtures
for j=1:NumSpectra
    %Use Monte Carlo method to determine abundances in each mixture
    PCTsMORB(j)=((HighAbunMORB-LowAbunMORB).*rand(1)+LowAbunMORB);
    %First half of spectra will be mare
    if j<NumSpectra/2
        PCTsMareMature(j)=((HighAbunReg-LowAbunReg).*rand(1)+LowAbunReg);
        PCTsMareImmature(j)=((HighAbunReg-LowAbunReg).*rand(1)+LowAbunReg);
        PCTsHighlandsImmature(j)=0;
        PCTsHighlandsMature(j)=0;
        %Second half of spectra will be highlands spectra
    else
        PCTsMareMature(j)=0;
        PCTsMareImmature(j)=0;
        PCTsHighlandsImmature(j)=((HighAbunReg-LowAbunReg).*rand(1)+LowAbunReg);
        PCTsHighlandsMature(j)=((HighAbunReg-LowAbunReg).*rand(1)+LowAbunReg);
    end
    PCTsPyroxene(j)=((HighAbunPyr-LowAbunPyr).*rand(1)+LowAbunPyr);
    %Select a random interpolated MORB spectrum to simulate  hydration
    RandWatAmt(j)=round(rand(1)*max(WatInterp));
    %Round to the nearest 1 ppm
    [M I]=min(abs(RandWatAmt(j)-WatInterp));
    %Actual water amount in mixture is morb abundance*interpolated MORB
    %spectrum used
    ActWatAmt(j)=WatInterp(I)*PCTsMORB(j);
    WtPercent(j,:)=[PCTsMORB(j),PCTsHighlandsMature(j),PCTsHighlandsImmature(j),PCTsMareMature(j),PCTsMareImmature(j),PCTsPyroxene(j)]; %Wt%
    if j<NumSpectra/2
        WtPercent(j,4)=1-sum(WtPercent(j,1:3))-sum(WtPercent(j,5:6)); % Fill in with mature mare so sum=1
    else
        WtPercent(j,2)=1-WtPercent(j,1)-sum(WtPercent(j,3:6)); % Fill in with mature highlands so sum=1
    end
    %Placing of SSAs, grain sizes (d) and densities (rho) into arrays
    SSAMORBChosen=TotalWatInterpESPAT(I,:).';
    AM=[SSAMORBChosen,SSAHighlandsMature,SSAHighlandsImmature,SSAMareMature,SSAMareImmature,SSAPyrox].';
    rho=1000*[densMORB,densReg,densReg,densReg,densReg,densPyrox].';
    d=[dMORB,dReg,dReg,dReg,dReg,dPyrox].';
    %Weighting of SSAs by cross-sectional area
    AV=AM(1:6,:)./(rho(1:6).*d(1:6))/(1/(mean(rho(1:6))*mean(d(1:6))));
    Numerator=WtPercent(j,:).'.*AM./(rho.*d);
    Denominator=WtPercent(j,:).'./(rho.*d);
    % Nonlinear mixing of SSAs to create intimate mixture via Hapke Radiative Transfer modeling
    SSAInt=sum(Numerator,1)/sum(Denominator); %SSA of intimate mixture
    % Need to convert to lidar geometry using SSA and scattering asymmetry
    % factor (p) of 1.5 (see manuscript for details on Hapke model)
    RInt=Hapke_Lidar_R_Function(SSAInt,1.5,WLS);
    SimulationLidarR(j,:)=RInt;
    SimulationSSA(j,:)=SSAInt;
end
%%
%Add random Gaussian noise to reflectance based on the instrument SNR
SNR=250;
NoiseVarIntimate=SimulationLidarR/SNR;
NoiseStd=sqrt(NoiseVarIntimate*1000)/1000;
SimulationLidarRNoisy=real(SimulationLidarR+NoiseStd.*randn(size(NoiseStd)));
%% Plot a subset of the mixtures in Reflectance
figure(2);
if NumSpectra<200
    plot(WLS,SimulationLidarRNoisy(:,:),'color',[0,0,0,0.5]);
else
    plot(WLS,SimulationLidarRNoisy(1:200,:),'color',[0,0,0,0.5]);
end
xlabel('Wavelength (\mum)','FontSize',16,'FontWeight','bold');
ylabel('Lidar Reflectance','FontSize',16,'FontWeight','bold');
ax=gca;
ax.LineWidth = 2;
ax.TickDir='in';
ax.FontName='arial';
ax.FontSize=14;
ax.FontWeight='bold';
box on
grid on
%% Using lidar retrieval to determine water abundance starts here
%Choose spectral (laser) wavelengths to measure with
LaserLams=1E-3*[1500,2650,2800,3100];
for k=1:numel(LaserLams)
    [val ind]=min(abs(LaserLams(k)-WLS));
    LLindices(k)=ind;
end
% Choose Constituents to Solve For
WaterPPM=[1522,762,176,22,0]; %Water and Hydroxyl from SIMS
MorbSelection1=1; %1522 ppm endmember
MorbSelection2=4; %22 ppm endmember
%Choose which MORB spectra to solve for
SSAMORB1=ConstituentSSAs(:,MorbSelection1);
SSAMORB2=ConstituentSSAs(:,MorbSelection2);
%Choose which other endmembers to solve for
AVolSolve=[SSAMORB1,SSAMORB2,SSAHighlandsImmature,SSAMareMature].';
rhoSolve=1000*[densMORB,densMORB,densReg,densReg].';
dSolve=[dMORB,dMORB,dReg,dReg].';
%Scale by grain size and density
AV=AVolSolve./(rhoSolve.*dSolve)/(1/(mean(rhoSolve)*mean(dSolve)));
%% Perform retrieval for each synthetic spectrum
for j=1:NumSpectra
    %Downselect to laser wavelengths
    WLSCut=LaserLams.';
    ObsSpectrum=SimulationLidarRNoisy(j,:);
    %Convert from noisy lidar reflectance to noisy SSA
    ObsSSA=Hapke_Lidar_SSA_function(ObsSpectrum,1.5,WLS);
    for i=1:numel(LLindices)
        %Select only measured and endmember SSAs at laser wavelengths
        AVol(:,i)=AV(:,LLindices(i));
        ObsIntimateSSACut(i)=ObsSSA(LLindices(i));
    end
    %Perform non-negative least squares algorithm using four endmember spectra
    %and observed SSA values of mixture
    MeasuredAbundances=lsqnonneg(AVol.',ObsIntimateSSACut.');
    MeasuredVolAmplitudes(j,:) = MeasuredAbundances;
    
    %Calculate total water amount using solved-for MORB endmember abundances at total
    %water amounts.
    AreaFactor=0.8428; %This factor is the ratio between the mean grain size and density difference
    %for the mixtures as generated (4 regolith, 1 Morb, 1 Pyroxene) and the retrieved
    %mixtures (2 regolith, 2 MORB.
    MeasuredWater(j,:)=MeasuredAbundances(1)/0.8428*WaterPPM(MorbSelection1)+MeasuredAbundances(2)/0.8428*WaterPPM(MorbSelection2); %retrieve two MORBs
    %Calculate total water error
    WatError(j,:) = MeasuredWater(j)-ActWatAmt(j);
end
%% Plot the results showing scatter around 1:1 retrieval and histogram of total water error
MareError=WatError(1:end/2);
HighlandsError=WatError(end/2+1:end);
figure(3);
subplot(2,1,1)
scatter(ActWatAmt(1:end/2),MeasuredWater(1:end/2),10,'filled','markerfacecolor','blue');
xlabel('Input Abundance (ppm)');
ylabel('Retr. Abundance (ppm)');
line([0,500],[0,500],'color','black','linestyle','--');
box on
grid on
title('Mare')
subplot(2,1,2)
histogram(MareError,50,'facecolor','blue');
xlabel('Total Water Error (ppm)');
ylabel('Frequency');
box on
grid on
title('Mare');

figure(4);
subplot(2,1,1)
scatter(ActWatAmt(end/2+1:end),MeasuredWater(end/2+1:end),10,'filled','markerfacecolor','red');
xlabel('Input Abundance (ppm)');
ylabel('Retr. Abundance (ppm)');
line([0,500],[0,500],'color','black','linestyle','--');
box on
grid on
title('Highlands');
subplot(2,1,2)
histogram(HighlandsError,50,'facecolor','red');
xlabel('Total Water Error (ppm)');
ylabel('Frequency');
box on
grid on
title('Highlands');