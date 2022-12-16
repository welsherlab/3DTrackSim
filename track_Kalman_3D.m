%General case Kalman track. to enter uniform background define r as empty matrix [] and enter only a single value in
%bofr
function [trkerr,len,tau,mleerr,photons,stg_out,part_out,kalm_out,N,sb,aoe]=track_Kalman_3D(D,s,bofr,r,ki,kiz,N,tau,varargin)
%Simulate 2D tracking with realstic stage response and uneven background
%intensity
%20323 KDW EDIT 200727 AJN
% DO NOT initialize StageCom X/Y or will convolve end of stg command/0
% array and stage will never move


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
%kalm_out - x, y, z position estimate relative to stage center
%N - max number of steps in trajectory
%sb - signal, and background photon generating CRs for entire trajectory.
%aoe - axis of escape coding. 1 for x, 2 for y, 3 for z, 4 for particle reached full
%duration. 

%signal background estimation flag
p = inputParser;
addParameter(p,'send',s,@isnumeric)%optional input to vary signal intensity over course of a traj using exponential decay. 
addParameter(p,'Dest',D,@isnumeric)%enables passing of misleading D for use in estimation while still generating trajectory with true diffusion coefficient. Defaults to truth. 
parse(p,varargin{:})
if p.Results.send>s
    error('Why would intensity increase over the course of a traj? Submit less stupid inputs.')
end
Dest=p.Results.Dest;

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

%define signal intensity for entirety of traj.
if s==p.Results.send
    s(1:N)=s;
else
    s(1:N)=s*((p.Results.send/s).^(1/double(N))).^(double(1:N));
end

%assorted initializations
stg_out=zeros(N,3);
photons=zeros(N,1);
kalm_out=zeros(N,3);
stageCom=zeros(N,3);
mnb_act=zeros(N,1);

%Integration
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
wz=(sigma_z).^2;%because z spread is ~2x that of in xy?

f = 70e3;     %TAG lens frequency
amp = 1e-6;   %TAG max scan amplitude (amplitude 1 micron total scan size 2)
binspertag=24; %this makes subtau/z binning 2x finer than the bins in xiaochen's existing sim which has 14 bins per xy scan pos. 
tp=1/f;
binsperz=round(binspertag*tau/tp);
tagtau=tau/binsperz;
b=zeros(1,binsperz);

%Initialize Kalman
x_k_1=[0;0];
z_k_1=0;
z_k=0;
sigma_xk_1=sigma_xy;
sigma_zk_1=sigma_z;
bspp=0;

%Size of KT
d=500e-9;

%Stage position
prev_stgxy=[0;0];
prev_stgz=0;

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
for k=1:int32(N)
    %First, get simulated observed number of photons based on the particle
    %and laser positions
    [xlaser,ylaser]=knightsTour(k,d);   %xy laser positions
    laserPosxy=[xlaser;ylaser];
    %pull tag phase for subbins
    laserPosz=tagLens((1+(k-1)*binsperz:k*binsperz),f,tagtau,amp);
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
    n=sum(nseq);
    photons(k)=n; %store total photons for out
    
    %Kalman comparison
    x_k = (sigma_xy^2*x_k_1+n*sigma_xk_1*laserPosxy)/(sigma_xy^2+n*sigma_xk_1);
    sigma_xk = sigma_xy^2*sigma_xk_1/(sigma_xy^2+n*sigma_xk_1);
    
    %Update Kalman variables for xy
    x_k_1=x_k;
    sigma_xk_1=sigma_xk+2*Dest*tau;
    
    %Do z Kalman if there is a photon along z
    %if there are no photons in a given bin, just add to bins since
    %previous photon
    if n~=0
        for subix=1:binsperz
            %increment bins since previous photon
            bspp=bspp+1;
            if nseq(subix)~=0 
                %kalman estimate along z
                z_k = (sigma_z^2*z_k_1+nseq(subix)*sigma_zk_1*laserPosz(subix))/(sigma_z^2+nseq(subix)*sigma_zk_1);
                sigma_zk = sigma_z^2*sigma_zk_1/(sigma_z^2+nseq(subix)*sigma_zk_1);
                
                %update step
                z_k_1=z_k;
                sigma_zk_1=sigma_zk+2*Dest*bspp*tagtau;
                
                %reset bins since previous photon
                bspp=0;
            end
        end
    else
        %if no photons were detected bins since previous photon increases
        %by total possible subbins
        bspp=bspp+binsperz;
    end
    
    kalm_out(k,:)=[x_k(1),x_k(2),z_k];
    
    %Integrate errors and Move "stage"
    if k>stepstowait
        %integrate error
        errint_x=errint_x+x_k(1);
        errint_y=errint_y+x_k(2);
        errint_z=errint_z+z_k;
        
        %very important axes must be rows time must be columns for fft conv
        stageCom(k,:)=[ki*errint_x, ki*errint_y, kiz*errint_z];

        stg_x=impulseResponse(stageCom(1:k,1),measuredImpulseResponse_glvo);
        stg_y=impulseResponse(stageCom(1:k,2),measuredImpulseResponse_glvo);
        if k<=neIR
            %original manual convolution
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
            stg_z=stg(lengthImpResp,:);          % Take just the most recent element
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
%% Calculate Outputs
len=k;
sb=[s',mnb_act];
% calculate err 
trkerr=[sqrt(sum((stg_out(:,1:2)-part_out(:,1:2)).^2,2)),abs(stg_out(:,3)-part_out(:,3)),sqrt(sum((stg_out(:,:)-part_out(:,:)).^2,2))];
mleerr=[sqrt(sum((kalm_out(1:len,1:2)+stg_out(1:len,1:2)-part_out(1:len,1:2)).^2,2)),abs(kalm_out(1:len,3)+stg_out(1:len,3)-part_out(1:len,3)),sqrt(sum((kalm_out(1:len,:)+stg_out(1:len,:)-part_out(1:len,:)).^2,2))];

 %% Make Plots
% Plot overall trajectory
figure(3)
clf
hold on
plot(part_out(:,1),'Color',[0.07,0.62,1],'LineWidth',2)
plot(stg_out(:,1),'Color',[0 0 1],'LineWidth',3)
% plot(kalm_out(:,1)+stg_out(:,1),'Color',[0.47 0.67 1])
plot(part_out(:,2),'Color',[0.85 0.33 0.10],'LineWidth',2)
plot(stg_out(:,2),'Color',[0.64 0.08 0.18],'LineWidth',3)
% plot(kalm_out(:,2)+stg_out(:,2),'Color',[0.99,0.77,0.26])
plot(part_out(:,3),'Color',[0.72 0.27 1],'LineWidth',2)
plot(stg_out(:,3),'Color',[0.27 0.04 0.42],'LineWidth',3)
% plot(kalm_out(:,3)+stg_out(:,3))
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

% Impulse Response Function
function out = impulseResponse(x,measuredImpulseResponse)
%Convolve with stage command data
temp = conv(x,measuredImpulseResponse);
%Crop data
out = temp(1:length(x));
end

function out=emRate(xy,z,cxy,cz,b,s,w,wz)
%x - position of particle
%c - position of beam
%b - background
%w - spatial covariance
%slower implementations matlab doesn't like inv construction.
% out=s*exp(-0.5*(xy-cxy)'*inv(w)*(xy-cxy))*exp(-0.5*(z-cz).^2/wz)+b;
% out=s*exp(-0.5*(xy-cxy)'/w*(xy-cxy))*exp(-0.5*(z-cz).^2/wz)+b;
out=s*exp(-0.5*(xy-cxy)'*(w\(xy-cxy)))*exp(-0.5*(z-cz).^2/wz)+b;
end


% function [y]=fconv(x, h)
% %FCONV Fast Convolution
% %Version 1.0
% %Original coded by: Stephen G. McGovern, 2003-2004.
% %Highly modified by: Stacey Niver 2022
% %converts stage command to columns
% Ly=length(x)+length(h)-1;  % 
% Ly2=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly
% X=fft(x, Ly2);             % Fast Fourier transform
% H=fft(h, Ly2);	           % Fast Fourier transform
% Y=X(:).*H(:);              % Multiplication in frequency space is convolution, : indexing ensures compatible orientation of vectors
% y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
% y=y(1:length(x));          % Take just the first N elements
% end

