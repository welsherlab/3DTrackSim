%General Case. to enter uniform background define r as empty matrix [] and enter only a single value in
%bofr. General case tracking with realistic stage response, uneven
%background intensities, and active sb estimation during tracking. Fixed kt
%scan pattern.
%bin time tau is limited to full TAG periods. 
%limiting bin time to 1.5 TAG periods would introduce asymmetry in pixel height
%sampling if one needed to predict bg seperately along xy and z. (which one
%might) such a scan of high and low pixels could be useful. 

%not true 1D, middle value from 2D bayes LUT must be used to scale LUT
%along z or code doesn't work. 
function [trkerr,len,tau,mleerr,photons,stg_out,part_out,posest_out,N,sbandests,aoe]=track_XYBayesZBayes(D,s,bofr,r,ki,kiz,N,ogtau,varargin)
%Simulate 2D tracking with realistic stage response and uneven background
%intensity
%20323 KDW EDIT 200805 AJN
%201117 AJN optimization of run time

%Inputs
%D - Diffusion coefficient in m^2/s
%s - particle brightness in cps. These are the number of counts expected
%when the particle is under constant illumination
%bofr - 1D vector 1 longer than length of r corresponding to background
%levels around r values. Selected based on first value radius is less
%than,final value in r is point at or above which bg is las value.
%r - 1D vector of all points at which bg changes.
%ki - integral feedback parameter along xy axes. 
%kiz - integral feedback parameter along z axis.
%N - initial proposed number of steps in trajectory
%tau - initial proposed bin time

%Outputs
%len - length of trajectory in steps
%tau - final trajectory bin time after adjustments
%trkerr - (xy, z, xyz) absolute difference between particle position and stage position
%mleerr - (xy, z, xyz) absolute difference between particle position and estimated
%position, data prior to 22/08/12 has only two columns (xy, z). 
%photons - observed photons for each bin
%stg_out - x, y, z stage position/laser scan center
%part_out - x, y, z particle positions
%posest_out - x, y, z position estimate relative to stage center
%N - max number of steps in trajectory
%sbandests - signal, background, signal estimate, and background estimates
%for entire trajectory. 
%aoe - axis of escape coding. 1 for x, 2 for y, 3 for z, 4 for particle reached full
%duration. 

%signal background estimation flag
p = inputParser;
addParameter(p,'sbest',true,@islogical)%optional input to toggle signal and bg estimation off. Defaults to true.
addParameter(p,'send',s,@isnumeric)%optional input to vary signal intensity over course of a traj using exponential decay. 
addParameter(p,'Dest',D,@isnumeric)%enables passing of misleading D for use in estimation while still generating trajectory with true diffusion coefficient. Defaults to truth. 
parse(p,varargin{:})
if p.Results.send>s
    error('Why would intensity increase over the course of a traj? Submit less stupid inputs.')
end
Dest=p.Results.Dest;

%TAG and EOD Scan Parameters.
d=500e-9;       %Size of KT
f = 70e3;       %TAG lens frequency
amp = 1e-6;     %TAG max scan amplitude (amplitude 1 micron total scan size 2)
binspertag=24;    %this makes subtau/z binning 2x finer than the bins in xiaochen's existing sim which has 14 bins per xy scan pos.
%bins per z recalculated later based on new bin time, this is bins per 20
%us

%check inputs
if (length(bofr)-1)~=length(r)
    error('input dimensions of r and background at given r are incompatible')
end
if length(r)>1
    if any(diff(r)<=0)
        error('radius of background flips must increase monotonically to avoid ambiguous definitions')
    end
end
if N>=2^31-1
    error('trajectory duration exceeds step numbers numerical precision, increase precision of N or shorten trajectory')
end

%if no radius for flip define arbitrary radius r and same background on
%both sides of flip
if isempty(r)
    bofr=[bofr,bofr];
    r=1;
end

%Round tau to nearest TAG period which introduces reproducible photon
%indexing/binning. 
tp=1/f;
% tau=tp*1.5*round(ogtau/tp/1.5);
tau=tp*1*round(ogtau/tp/1);
if tau<=20e-6 %computation cannot compile/run in less than 20 so min bin size is 20
    tau=tp*1*ceil(ogtau/tp/1);
end
gsmax=41;
%Define Preliminary grid size
if Dest~=0
    gs=ceil(1e-6/1.5/sqrt(2*Dest*tau));
else
    gs=gsmax;
end
%ensure laser scan positions fall exactly onto grid by fixing grid size to
%nearest multiple of 5
if rem(gs-5,4)==0
elseif rem(gs-5,4)==1
    gs=gs-1;
elseif rem(gs-5,4)==2
    gs=gs+2;
elseif rem(gs-5,4)==3
    gs=gs+1;
end
%limit computation duration to accurate FPGA implementation max grid size by recalculating tau to the nearest us
if gs>gsmax
%     warning('grid too large for FPGA implementation at current tau')
    %recalculate tau rounded to the nearest ms
    tau=round((1e-6/1.5/gsmax)^2/2/Dest,6);
    %Define finalized modified grid size if necessary
    gs=ceil(1e-6/1.5/sqrt(2*Dest*tau));
    %ensure laser scan positions fall exactly onto grid
    if rem(gs-5,4)==0
    elseif rem(gs-5,4)==1
        gs=gs-1;
    elseif rem(gs-5,4)==2
        gs=gs+2;
    elseif rem(gs-5,4)==3
        gs=gs+1;
    end
else
    %do not change tau from the rounded to nearest TAG period value
end
if gs>gsmax
    error('correction by changing tau failed so the round vs ceil thing did matter rethink this implementation')
end
if gs<5
    error('particle motion is too fast to be tracked without expanding scan area')
end
mdptindx=ceil(gs/2);
%finalize tau by rounding to nearest TAG Period
% tau=tp*1.5*round(tau/tp/1.5);
tau=tp*1*round(tau/tp/1);
if tau~=ogtau
    %recalculate N to preserve intended trajectory duration.
    N=int32(ceil(N*ogtau/tau));
    if tau<ogtau
%         disp(['Tau has been decreased to ' num2str(tau) ' N has been increased to ' num2str(N) ' to preserve intended trajectory duration'])
    else
%         disp(['Tau has been increased to ' num2str(tau) ' N has been decreased to ' num2str(N) ' to preserve intended trajectory duration'])
    end
end

%after N has been recaulculated define signal intensity for traj.
if s==p.Results.send
    s(1:N)=s;
else
    s(1:N)=s*((p.Results.send/s).^(1/double(N))).^(double(1:N));
end

binsperz=round(binspertag*tau/tp);
%after final tau established calculate number of bins per EOD bin.
tagtau=tau/binsperz;
b=zeros(1,binsperz); %initialize b values for photon generation.

%Generate XY and Z coordinates, granularity along z is 2 x gs of xy to make
%voxels square making D universal instead of axially dependent.
coord=linspace(-0.5e-6,0.5e-6,gs);
[x2D,y2D]=meshgrid(coord);
z1D=linspace(-amp,amp,2*gs-1);%z axis scan is 2x xy so equivalent grid spacing will need 2x pts. 
% z1D=linspace(-amp,amp,501);%give nice fine resolution for dev purposes

%Quality control bins per z. sampling must have at least 8 pts/sine wave
%and meet or exceed bins along z axis.
if binsperz<tau/tp*8
    error('TAG sampling is inadequate to reproduce properties of sine wave, TAU must be subdivided into a minimum of 8 points per period to ensure adequate sampling, increase points per subbin')
end

%Initialize priors
Pr_xy=single(ones(size(x2D)));       %2D
Pr_z=single(ones(size(z1D)));

%Normalize Priors
Pr_xy=Pr_xy/sum(Pr_xy(:));       %2D
Pr_z=Pr_z/sum(Pr_z);

%Log Prior
lPr_xy=log(Pr_xy);            %2D
lPr_z=log(Pr_z);

%Calculate the diffusion convolution kernels along xy and z seperately
%this simulation has fixed them such that they are the same in the
%definition of the coordinate system along z however the odds of me
%changing that at some point are high so for reasons of robustness the
%implementation here stands.
if Dest~=0
    %XY diffusion Kernel Generation
    difffull=exp(-coord.^2/4/Dest/tau);
    c=find(difffull>=.1*max(difffull));
    %check symmetry of coordinates
    mid=ceil(length(difffull)/2);
    if mid-c(1)==c(end)-mid
        diffkernxy=difffull(c);
    else
        warning('Diffusion kernel is asymmetric along xy')
        diffkernxy=difffull(c);
    end
    %Z diffusion Kernel Generation
    difffull=exp(-z1D.^2/4/Dest/tau);
    c=find(difffull>=.1*max(difffull));
    %check symmetry of coordinates
    mid=ceil(length(difffull)/2);
    if mid-c(1)==c(end)-mid
        diffkernz=difffull(c);
    else
        warning('Diffusion kernel is asymmetric along z')
        diffkernz=difffull(c);
    end
    %quality control diffusion coefficients
    if length(diffkernxy)~=3%||length(diffkernz)~=3
        error('diffusion kernel size is larger than 3 pts')
    end
    
    if abs(diffkernxy(1)-0.3247)/0.3247>0.5%||abs(diffkernz(1)-0.3247)/0.3247>0.5
        error('Tau modifications have led to deviation from targeted diffusion kernel by more than 50%, actual diffusion kernel value falls outside of range from 0.1623 to 0.4870')
    end
%     disp(['zaxis diffusion kernel is ' num2str(length(diffkernz)) ' units long'])
    clear c mid
end

%assorted output initializations
stg_out=zeros(N,3);
photons=zeros(N,1);
posest_out=zeros(N,3);
stageCom=zeros(N,3);
sb_est=zeros(N,2);
mnb_act=zeros(N,1);

%Integration time
errint_x=0;
errint_y=0;
errint_z=0;

%Define steps to wait before turning on tracking
stepstowait=0;

%Generate Trajectory
%Generate steps in x (dx) and y (dy)
dx = sqrt(2*D*tau) * randn(N,1)';
dy = sqrt(2*D*tau) * randn(N,1)';
dz = sqrt(2*D*tau) * randn(N,1)';

%Generate XY trajectory
x_part = cumsum(dx)';
y_part = cumsum(dy)';
z_part = cumsum(dz)';
part_out=[x_part,y_part,z_part];

%PSF size
sigma_xy=128e-9;
sigma_z=718e-9;
%PSF Measurements empirical Chen Zhang https://doi.org/10.3390/e23050498
w=sigma_xy^2*eye(2,2);
wz=(sigma_z).^2;

%precalculate partial gamma to generate LUT
partialgamma_xy=single(zeros(length(coord),length(coord),25));
for i=1:25
    [xlaser,ylaser]=knightsTour(i,d);
    partialgamma_xy(:,:,i)=exp(-((x2D-xlaser).^2./2./sigma_xy./sigma_xy)-((y2D-ylaser).^2./2./sigma_xy./sigma_xy));
end

% %partial gamma along z LUT
% %calculate p(n) photons for all grid positions assuming laser exists at
% %all grid positions
% partialgamma_z=zeros(binsperz,length(z1D),2);
% for k=1:2
% laserPosz=tagLens((1+(k-1)*binsperz:k*binsperz),f,tagtau,amp);
% partialgamma_z(:,:,k)=exp(-((z1D-laserPosz').^2./2./sigma_z./sigma_z));
% end

%Stage position
prev_stgxy=[0;0];
prev_stgz=0;

% Initialize photon characterization for signal and background estimation
pixA=ismember(mod((1:N)-1,25)+1,23);
pixF=ismember(mod((1:N)-1,25)+1,[1,3,5,7]);
testimate=15e-3; %rolling average of 15 ms worth of data goes into s and b estimation

binsestimate=int32(testimate/tau); %number of bins that are equal to 15 ms of data
numA=sum(pixA(1:binsestimate));
numF=sum(pixF(1:binsestimate));
pixA=[];
pixF=[];

% original cal dist weight (σXY = 0.038 µm, σZ = 0.15 µm).
%"Z:\Stacey Z\MATLAB\TrackingSimulations\220519_3DTrack\220726_AFweight_StoOCRconversion"
% weight=[0.5905,1;1.9853e-7,1];
%Z:\Stacey Z\MATLAB\TrackingSimulations\220519_3DTrack\221207_AFweight_UsingExperimentallyObservedSigmasFromEvenTrk
weight =[0.44750577	1
1.0185348e-05	1];

% Load impulse response data
measuredImpulseResponse_stg = dlmread('Impulse_Response.txt');
load('galvoimpulseresponse.mat','y')
measuredImpulseResponse_glvo = y;
if tau==20e-6
    %do nothing undecimated Impulse Response Function is Entirely adequate
elseif rem(tau/20e-6,1) == 0
    %if integer relationship decimate measured impulse response to
    %appropriate sampling rate
    measuredImpulseResponse_stg = decimate(measuredImpulseResponse_stg,tau/20e-6);
    measuredImpulseResponse_glvo = decimate(measuredImpulseResponse_glvo,tau/20e-6);
else
    %interpolate to nearest us sampling then decimate
%     warning('tau was rounded to nearest us for resampling stage response')
    measuredImpulseResponse_stg = interp(measuredImpulseResponse_stg,20);
    measuredImpulseResponse_stg = decimate(measuredImpulseResponse_stg,double(int32(round(tau,6)/1e-6)));
    measuredImpulseResponse_glvo=interp(measuredImpulseResponse_glvo,20);
    measuredImpulseResponse_glvo=decimate(measuredImpulseResponse_glvo,double(int32(round(tau,6)/1e-6)));
end
neIR=length(measuredImpulseResponse_stg);
if neIR<1000
    warning('Fourier Convolution May no Longer Be Temporally Efficient/Faster than traditional conv')
end

%precalculate fft of measured Impulse Response
Ly=length(measuredImpulseResponse_stg)*2-1;  %
ftzerpad=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly
fftExpImpResp=fft(measuredImpulseResponse_stg, ftzerpad);             % Fast Fourier transform
lengthImpResp=length(measuredImpulseResponse_stg);

%Now simulate tracking
for k=int32(1:N)
    %First, get simulated observed number of photons based on the particle
    %and laser positions
    [xlaser,ylaser]=knightsTour(k,d);   %xy laser positions
    laserPosxy=[xlaser;ylaser];
    %pull tag phase for subbins
    laserPosz=tagLens((1+(k-1)*binsperz:k*binsperz),f,tagtau,amp);
    if k==1 %on first loop check relative spacing.
        edges=linspace(-1e-6,1e-6,3);%only need two bins high and low n, less bins less math. 
        if mean(diff(unique(laserPosz)))>=mean(diff(edges))
            warning('photons are being generated at a lower z resolution than they are being assigned to positions for bayesian estimation')
        end
        [nbinperk,~,bin]=histcounts(laserPosz,edges);
        bin=bin';%indexes must be column for accumarray.
        nbinperk=nbinperk'; %pretty sure I need this to be a column also for consistency
        laserPosz_binned=accumarray(bin,laserPosz)./nbinperk;
        partialgamma_z=exp(-(z1D-laserPosz_binned).^2./2./sigma_z/sigma_z);
    end
    %radius of laser position from initial stage position 0,0 to determine
    %background intensity
    laser_r=sqrt(sum((prev_stgxy+laserPosxy).^2)+laserPosz.^2);
    %particlePosxy for given step in correct matrix format
    particlePosxy=[x_part(k);y_part(k)];
    
    %Determine actual bg at given laser positions
    for zscn=1:binsperz
        if laser_r(zscn)>r(end)
            %if laser radius is outside last point of bg increase, set to value
            %at furthest radius
            b(zscn)=bofr(end);
        else
            %otherwise set background to appropriate level
            b(zscn)=bofr(find(laser_r(zscn)<r,1,'first'));
        end
    end
    mnb_act(k)=mean(b);
    
    %calculate observed number of photons
    nseq=poissrnd(tagtau.*emRate((particlePosxy-prev_stgxy),(z_part(k)-prev_stgz),laserPosxy,laserPosz,b,s(k),w,wz)); %number of photonsd
    nseq=nseq';
    n=sum(nseq);
    photons(k)=n; %store total photons for out
    ktindx=mod(k-1,25)+1;
    
    %if photons belong to pixels of interest append them to photon string
    if ismember(ktindx,23)
        pixA=[pixA,photons(k)];
    elseif ismember(ktindx,[1,3,5,7])
        pixF=[pixF,photons(k)];
    end
    
    %estimate signal and bg levels based on observed photons
    if k<=binsestimate
        %some given estimate. Can be accurate or could be not. We will make
        %it accurate. But at less than 15 ms of data that's hard to deal
        %with
        [calc_sb]=[s(k);mnb_act(k)];
    else
        %above 15 ms of available data estimate using only the most recent
        %15 ms worth of data.
        %calculate index value for end of given photon type string
        %then select most recent desired number of photons for a given sb
        %estimation time window calculated early
        A=length(pixA);
        photonsA=pixA(A-numA+1:A);
        
        F=length(pixF);
        photonsF=pixF(F-numF+1:F);
        
        [calc_sb]=weight\[sum(photonsA)/(length(photonsA)*tau);sum(photonsF)/(length(photonsF)*tau)];
    end
    if p.Results.sbest==false
        %give perfect info
        [calc_sb]=[s(k);mnb_act(k)];
    end
    %if signal estimate is negative or 0 crash
        if calc_sb(1)<=0
            warning('negative or 0 estimated signal')
            break
        end
    
        %if bg is negative replace it with 0
        if calc_sb(2)<0
            calc_sb(2)=0;
        end
    
        %store signal/background estimates for a given step along with the
        %actual bg level
        sb_est(k,:)=calc_sb';
     
    %% use 2D Bayesian to Estimate Particle Position Along X and Y
    %Construct likelihood of detection n photons. This will be a 2D
    %code run speed is increased by combining arguments/doing single exp
    %gamma2D=calc_sb(1)*exp(-(x2D-xlaser).^2./2./sigma./sigma).*exp(-(y2D-ylaser).^2./2./sigma./sigma)+calc_sb(2);
    %code run speed is further improved by recalculating a partial gamma
    %LUT for each possible laser scan position in a 25 pt kt. Does not
    %apply to non kt scan pattern.
    %gamma2D=calc_sb(1)*exp(-((x2D-xlaser).^2./2./sigma./sigma)-((y2D-ylaser).^2./2./sigma./sigma))+calc_sb(2);
    gammaxy=calc_sb(1)*partialgamma_xy(:,:,ktindx)+calc_sb(2);
    
    %The non-log likelihood can lead to infs at high n. Try log
    %P_n_2D=((gamma2D*tau).^n).*exp(-gamma2D.*tau)./factorial(n);
    log_n_factorial=stirling(n);
    lP_n_xy=n*log(gammaxy*tau)-gammaxy*tau-log_n_factorial;
    
    %Calculate log posterior
    lPoxy=lPr_xy+lP_n_xy;
    
    %Now normalize log posterior using log-sum-exp method
    lPoxy=lPoxy-max(max(lPoxy));
    lPoxy=lPoxy-log(sum(sum(exp(lPoxy))));
    
    %Account for diffusion via convolution. If D=0, then posterior becomes
    %prior
    if Dest~=0
        %old convolution uses full kernel is slow
        %lPr2D=log(conv2(exp(lPo2D),diff(481:521,481:521),'same'));
        %new convolution takes advantage of symmetry and convolves along
        %columns then along rows
        lPr_xy=log(conv2(diffkernxy,diffkernxy,exp(lPoxy),'same'));
        lPr_xy(isinf(lPr_xy))=-745;
    else
        lPr_xy=lPoxy;
    end
    
    %% use Bayesian to Estimate Particle Position Along Z
    % bin nseq along laser scan using indeces from earlier.
    binnedn=accumarray(bin,nseq);
    %pull the photons from the most recent n bins
    log_nz_factorial=stirling(binnedn);
    %gamma z taking into account kt scan position.
%     gammaz=calc_sb(1)*exp(-((xlaser).^2./2./sigma_xy./sigma_xy)-((ylaser).^2./2./sigma_xy./sigma_xy))*exp(-(z1D-laserPosz_binned).^2./2./sigma_z/sigma_z)+calc_sb(2);
    % for LUT version will need to pull middle value from gammaxypartial
    % for given kt index and laser pos z LUT. LUT version saves only like
    % 0.3 of second but I already wrote it so meh. 
    gammaz=calc_sb(1)*partialgamma_xy(mdptindx,mdptindx,ktindx)*partialgamma_z+calc_sb(2);
    
    
    %calculate log(P(n)) and posterior
    lP_n_z=binnedn.*log(tagtau.*nbinperk.*gammaz)-gammaz.*tagtau.*nbinperk-log_nz_factorial;
    lPoz=lPr_z+sum(lP_n_z);
    
    %normalize using log sum exp
    lPoz=lPoz-max(lPoz);
    lPoz=lPoz-log(sum(exp(lPoz)));
    
    %Account for diffusion via convolution. If D=0, then posterior becomes
    %prior
    if Dest~=0
        %old convolution uses full kernel is slow
        %lPr2D=log(conv2(exp(lPo2D),diff(481:521,481:521),'same'));
        %new convolution takes advantage of symmetry and convolves along
        %columns then along rows
        lPr_z=log(conv(exp(lPoz),diffkernz,'same'));
        lPr_z(isinf(lPr_z))=-745;
    else
        lPr_z=lPoz;
    end

%     %reset binned z arrivals
%     binnedsubn=zeros(length(z1D_scnsubset),1);
    
    %% Find maximum likelihood estimate of position
    temp=x2D(lPoxy==max(max(lPoxy)));
    posest_out(k,1)=mean(temp);
    temp=y2D(lPoxy==max(max(lPoxy)));
    posest_out(k,2)=mean(temp);
    temp=z1D(lPoz==max(lPoz));
    posest_out(k,3)=mean(temp);
    
    if isnan(posest_out(k,1))||isnan(posest_out(k,2))
        warning('return of the NaNs - hunt until dead')
        break
    end
    
    %Integrate errors and Move "stage"
    if k>stepstowait
        %integrate error
        errint_x=errint_x+posest_out(k,1);
        errint_y=errint_y+posest_out(k,2);
        errint_z=errint_z+posest_out(k,3);
        
        %very important axes must be rows time must be columns for fft conv
        stageCom(k,:)=[ki*errint_x, ki*errint_y, kiz*errint_z];
        
        stg_x=impulseResponse(stageCom(1:k,1),measuredImpulseResponse_glvo);
        stg_y=impulseResponse(stageCom(1:k,2),measuredImpulseResponse_glvo);
        if k<=neIR
            stg_z=impulseResponse(stageCom(1:k,3),measuredImpulseResponse_stg);
            stg_z=stg_z(k);
        else
            %Fourier Based Convolution
            %between x (the stage command) and h the already fourier transpormed
            %measured impulse response using zp the calculated zero padding from
            %initial transformation. Critical that x be stage commands as Columns!
            fftStgCom=fft(stageCom(k-neIR+1:k,3), ftzerpad);             % Fast Fourier transform
            Y=fftStgCom.*fftExpImpResp;                    % Multiplication in frequency space is convolution
            stg=real(ifft(Y, ftzerpad));      % Inverse fast Fourier transform
            stg_z=stg(lengthImpResp,:);          % Take just the first N elements        end
        end
        stg=[stg_x(k) stg_y(k) stg_z];
        stg_out(k,:)=stg;
        prev_stgxy=[stg(1);stg(2)];
        prev_stgz=stg(3);
    end
    aoe=4;
    %if particle escaped box end trajectory
    if abs(prev_stgxy(1)-x_part(k))>(d+0.5e-6)
        aoe=1;
%         disp('particle escaped box along x')
        break
    end
    
    if abs(prev_stgxy(2)-y_part(k))>(d+0.5e-6)
        aoe=2;
%         disp('particle escaped box along y')
        break
    end
    
    if abs(prev_stgz-z_part(k))>(amp+.5e-6)
        aoe=3;
%         disp('particle escaped box along z')
        break
    end
end
%% calculate/define outputs
len=k;
sbandests=[s',mnb_act, sb_est];
% calculate err
trkerr=[sqrt(sum((stg_out(:,1:2)-part_out(:,1:2)).^2,2)),abs(stg_out(:,3)-part_out(:,3)),sqrt(sum((stg_out(:,:)-part_out(:,:)).^2,2))];
mleerr=[sqrt(sum((posest_out(1:len,1:2)+stg_out(1:len,1:2)-part_out(1:len,1:2)).^2,2)),abs(posest_out(1:len,3)+stg_out(1:len,3)-part_out(1:len,3)),sqrt(sum((posest_out(1:len,:)+stg_out(1:len,:)-part_out(1:len,:)).^2,2))];
%% Make Plots
% %
figure(3)
clf
hold on
plot(part_out(:,1),'Color',[0.07,0.62,1],'LineWidth',2)
plot(stg_out(:,1),'Color',[0 0 1],'LineWidth',3)
% plot(posest_out(:,1)+stg_out(:,1),'Color',[0.47 0.67 1])
plot(part_out(:,2),'Color',[0.85 0.33 0.10],'LineWidth',2)
plot(stg_out(:,2),'Color',[0.64 0.08 0.18],'LineWidth',3)
% plot(posest_out(:,2)+stg_out(:,2),'Color',[0.99,0.77,0.26])
plot(part_out(:,3),'Color',[0.72 0.27 1],'LineWidth',2)
plot(stg_out(:,3),'Color',[0.27 0.04 0.42],'LineWidth',3)
% plot(posest_out(:,3)+stg_out(:,3))
hold off
% legend('laserradius', 'particle x','stage x','kalman x','particle y','stage y','kalman y','particle z','stage z', 'kalman z')
legend( 'particle x','stage x','particle y','stage y','particle z','stage z')
xlabel('Time/tau')

figure(2)
clf
histogram(trkerr(:,1)*1e6)
hold on
histogram(trkerr(:,2)*1e6)
legend('xy','z')
xlabel('error in stage position (um)')

figure(4)
clf
histogram(mleerr(:,1)*1e6)
hold on
histogram(mleerr(:,2)*1e6)
hold off
legend('xy','z')
xlabel('error in position estimate (um)')
end

%Knight's tour coords
function [xOut,yOut]=knightsTour(index,d)
x=d*[-1;-0.5;-1;0;1;0.5;1;0;0.5;1;0;-1;-0.5;0.5;1;0.5;-0.5;-1;0;1;0.5;-0.5;0;-0.5;-1];
y=d*[1;0;-1;-0.5;-1;0;1;0.5;-0.5;0.5;1;0.5;-0.5;-1;0;1;0.5;-0.5;-1;-0.5;0.5;1;0;-1;0];
xOut=x(mod(index-1,length(x))+1);
yOut=y(mod(index-1,length(x))+1);
end

% TAG lens coordinates
function [zOut]=tagLens(index,f,tau,amp)
phi = cumsum(2*pi*f*tau*ones(1,round(1/(f*tau))));
z = amp * sin([0 phi(1:end-1)]);
zOut = z(mod(index-1,length(z)) + 1);
end

%Stirling approximation for log(n!)
function out=stirling(n)
out=zeros(size(n));
out(n<10)=log(factorial(n(n<10)));
out(n>=10)=n(n>=10).*log(n(n>=10))-n(n>=10)+1;
end

% Impulse Response Function
function out = impulseResponse(x,measuredImpulseResponse)
%Convolve with stage command data
temp = conv(x,measuredImpulseResponse);
%Crop data
out = temp(1:length(x));
end

function out = emRate(xy,z,cxy,cz,b,s,w,wz)
%x - position of particle
%c - position of beam
%b - background
%s - signal
%w - spatial covariance
out=s*exp(-0.5*(xy-cxy)'*(w\(xy-cxy)))*exp(-0.5*(z-cz).^2/wz)+b;
end
