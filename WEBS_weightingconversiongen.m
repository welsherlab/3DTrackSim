s=logspace(3,6.4771,30)./.09; %actual signal from experimental input calibration curve slope
b=0;
trials=50;
N=5e3;
tbin=20e-6;

% %original values
% sigmaxy = 0.38e-7;
% sigmaz = 1.5e-7;

% % observed COBWEBS-KALMAN values
% sigmaxy = 0.068e-6;
% sigmaz = 0.16e-6;
% %set filename
% fname='221207_newpixelweightsfor3Dtracks_CBWBSKAL_rawdata.mat';
% fnameb='221207_newpixelweightsfor3Dtracks_CBWBSKAL_weights.mat';

% observed COBWEBS-COBWEBS values
sigmaxy = 0.085e-6;
sigmaz = 0.12e-6;
%set filename
fname='221207_newpixelweightsfor3Dtracks_CBWBSCBWBS_rawdata.mat';
fnameb='221207_newpixelweightsfor3Dtracks_CBWBSCBWBS_weights.mat';

%even if we can get better bg estimations later on we don't particularly care bc beyond 100 ms it's definitely not useful? Seems fair to me. 
%Note N MUST be an even multiple of the number of positions in your scan
%pattern or the indexing will be incorrect and problems will ensue. In case
%of kt must be even multiple of 25. 
%Code runs super fast. consider doing 1000 trials and 10000 duration just
%to have a far out for suresies stable comparison point. 



if exist(fname,'file') ~= 2 %data does not already exist and must be generated
%Generate photons counts for particle randomly offset from stage center
%consistend with observed distribution from center of stage in simulated
%tracking from mean standard deviation of 200620 scattering data set and
%200726 in cell high data set with assumed mean of 0 which was close to
%approximation. 
%initialize variables
photons=zeros(length(s),trials,N,'single');
xc=zeros(1,N,'single');
yc=zeros(1,N,'single');

% data is stored as (signal level, background level, trial index, tracking
% index/time bin)
for i=1:length(s)
    for t=1:trials
        if i==1&&t==1
            %if is first time through loop define XC and YC since laser
            %positions are consistent for all trials
            [photons(i,t,:),xc,yc]=photongen3D(s(i),b,N,tbin,sigmaxy,sigmaz);
        else
            %subsequent trials do not save/export laser positions
            [photons(i,t,:),~,~]=photongen3D(s(i),b,N,tbin,sigmaxy,sigmaz);
        end
        
    end
    disp(['Signal set ' num2str(i) '/' num2str(length(s))])
end
save(fname)
else
    load(fname)
end

%now have data generate new photon sb conversion weights.
%indexing
%Knights tour divide photon counts into equivalent pixel pools lettered
%from the middle out
pixA=ismember(mod((1:N)-1,25)+1,23);
pixB=ismember(mod((1:N)-1,25)+1,[2,4,6,8]);
pixC=ismember(mod((1:N)-1,25)+1,[9,13,17,21]);
pixD=ismember(mod((1:N)-1,25)+1,[11,15,19,25]);
pixE=ismember(mod((1:N)-1,25)+1,[10,12,14,16,18,20,22,24]);
pixF=ismember(mod((1:N)-1,25)+1,[1,3,5,7]);

%Generate Index values corresponding to each type of pixels (1 per scan (A), 4
%per scan (B,C,D,F), and 8 per scan(E))
indx1=1:length(photons(pixA));
indx4=1:length(photons(pixF));
indx8=1:length(photons(pixE));
indx25=1:N;

% data is initially stored as (signal level, trial index, tracking
% index/time bin)
% permute photon data to put tracking index/bin number first so that matlab
% can actually do math with it
photons=permute(photons,[3,1,2]);

%flip orientation of indx values to vertical
i=size(indx1);
if i(1)==1
    indx1=indx1';
    indx4=indx4';
    indx8=indx8';
    indx25=indx25';
end

% DETERMINE OBSERVED COUNT RATES: for all 6 equivalent pixel types alone
%indx*tbin is time of observation over selected photons. each photons value
%is a single bin hence linear increase of indexing. OCR=counts/timespent
%observing a given pixel. cps just like s and b
obs_cr_A=cumsum(photons(pixA,:,:),1)./(tbin*indx1);
obs_cr_B=cumsum(photons(pixB,:,:),1)./(tbin*indx4);
obs_cr_C=cumsum(photons(pixC,:,:),1)./(tbin*indx4);
obs_cr_D=cumsum(photons(pixD,:,:),1)./(tbin*indx4);
obs_cr_E=cumsum(photons(pixE,:,:),1)./(tbin*indx8);
obs_cr_F=cumsum(photons(pixF,:,:),1)./(tbin*indx4);
obs_cr_all=cumsum(photons(:,:,:),1)./(tbin*indx25);

weight=[
    mean(obs_cr_A(end,:,:)./s,'all'),1;
    mean(obs_cr_F(end,:,:)./s,'all'),1;
    mean(obs_cr_all(end,:,:)./s,'all'), 1];

disp('weighting factors calculated')
save(fnameb,'weight')

apw=[mean(obs_cr_A(end,:,:)./s,'all')
    mean(obs_cr_B(end,:,:)./s,'all')
    mean(obs_cr_C(end,:,:)./s,'all')
    mean(obs_cr_D(end,:,:)./s,'all')
    mean(obs_cr_E(end,:,:)./s,'all')
    mean(obs_cr_F(end,:,:)./s,'all')];
pix2D=[
apw(6)  apw(5)  apw(4)  apw(5)  apw(6)
apw(5)  apw(3)  apw(2)  apw(3)  apw(5)
apw(4)  apw(2)  apw(1)  apw(2)  apw(4)
apw(5)  apw(3)  apw(2)  apw(3)  apw(5)
apw(6)  apw(5)  apw(4)  apw(5)  apw(6)];

A=figure;
imagesc(pix2D)
axis('off')
axis square
colorbar
saveas(A,'OverallSWeightPerPix.fig')
ax=gca;
ax.Children.AlphaData=[1 0 0 0 1;0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 1 0 0 0 1];

%Make plots of interest
%weighting factor histograms
A=figure (1);
histogram(obs_cr_A(end,:,:)./s)
title('A pixels-Signal Weighting')
xlabel('observed count rate/input count rate in 0 bg')
ylabel('# observations')
xline(weight(1,1));
annotation('textbox','String',['Weight = ' num2str(weight(1,1))],'FitBoxToText','on');
saveas(A,'s_weighting_A.fig')

A=figure (2);
histogram(obs_cr_F(end,:,:)./s)
title('F pixels-Signal Weighting')
xlabel('observed count rate/input count rate in 0 bg')
ylabel('# observations')
xline(weight(2,1));
annotation('textbox','String',['Weight = ' num2str(weight(2,1))],'FitBoxToText','on');
saveas(A,'s_weighting_F.fig')

%input s to OCR conversion for 3D plot gen
A=figure(3);
histogram(obs_cr_all(end,:,:)./s)
title('OCR/Sig full scan area')
xlabel('observed count rate/input count rate in 0 bg')
ylabel('# observations')
xline(weight(2,1));
annotation('textbox','String',['Weight = ' num2str(weight(3,1))],'FitBoxToText','on');
saveas(A,'s_weighting_fullscan.fig')

