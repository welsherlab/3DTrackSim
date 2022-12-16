%General case Kalman track. to enter uniform background define r as empty matrix [] and enter only a single value in
%bofr
function [photons,xc,yc]=photongen3D(s,b,N,tau,sigmaxy,sigmaz)
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
%mu_x - initial X particle position
%mu_y - initial Y particle position
%ki - integration constant for integral
%kp - integration constant for proportional
%N - number of steps in trajectory

%Outputs
%len - length of trajectory in steps
%trkerr - absolute difference between particle position and stage position
%photons - observed photons for each bin
%stg_out - x,y,z stage position
%part_out - x,y,z particle positions
%kalmx - x position estimate
%kalmy - y position estimate
%par_rad - particle radius

if N>=2^31-1
    error('trajectory duration exceeds step numbers numerical precision, increase precision of N or shorten trajectory')
end

%assorted initializations
photons=zeros(N,1);
xc=zeros(N,1);
yc=zeros(N,1);

%Generate particle positions Displacements from mean over displacement
x_part = random ('Normal',0,sigmaxy,[1,N]);
y_part = random ('Normal',0,sigmaxy,[1,N]);
z_part = random ('Normal',0,sigmaz,[1,N]);

%PSF size
sigma_xy=128e-9;
sigma_z=718e-9;
%PSF Measurements empirical Chen Zhang https://doi.org/10.3390/e23050498
w=sigma_xy^2*eye(2,2);
wz=(sigma_z).^2;%because z spread is ~2x that of in xy?

f = 70e3;     %TAG lens frequency
amp = 1e-6;   %TAG max scan amplitude (amplitude 1 micron total scan size 2)
binsperz=30; %this makes subtau/z binning 2x finer than the bins in xiaochen's existing sim which has 14 bins per xy scan pos.
tagtau=tau/binsperz;

%Size of KT
d=500e-9;

%Now simulate tracking
for k=1:int32(N)
    %First, get simulated observed number of photons based on the particle
    %and laser positions
    [xc(k),yc(k)]=knightsTour(k,d);   %xy laser positions
    laserPosxy=[xc(k);yc(k)];
    %pull tag phase for subbins
    laserPosz=tagLens((1+(k-1)*binsperz:k*binsperz),f,tagtau,amp);
    %particlePosxy for given step in correct matrix format
    particlePosxy=[x_part(k);y_part(k)];
    
    %calculate observed number of photons
    nseq=poissrnd(tagtau.*emRate((particlePosxy),(z_part(k)),laserPosxy,laserPosz,b,s,w,wz)); %number of photonsd
    n=sum(nseq);
    photons(k)=n; %store total photons for out
end
end

%Knight's tour coords
function [xOut,yOut]=knightsTour(index,d)
x=d*[ -1;-0.5;-1;0;1;0.5;1;0;0.5;1;0;-1;-0.5;0.5;1;0.5;-0.5;-1;0;1;0.5;-0.5;0;-0.5;-1];
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

