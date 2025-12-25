function [x, z, migRF] = fkmig_monostat(RF, fs, tstart, pitch, vsound, upsamp)
%FKMIG_MONOSTAT Uses Stolt's f-k Migration on Monostatic Synthetic Aperture
% [x, z, migsig] = fkmig_monostat(signal, fs, tstart, pitch)
% Inputs:
%   RF -- Signal array (row dimension -- time; column dimension -- lateral)
%   fs -- Sampling Rate Along Time Dimension [Hz]
%   tstart -- Time for First Sample in Array [s]
%   pitch -- Lateral Spacing Between Array Elements [m]
%   vsound -- Speed of Sound [m/s]
%   upsamp -- Vector of Upsampling in Lateral(1) and Axial(2) Dimensions
% Outputs:
%   x -- Lateral Dimension of Migrated RF Image
%   z -- Axial Dimension of Migrated RF Image
%   migRF -- Migrated RF Image

% Assembling Inputs
[nt0, nx0] = size(RF); % Number of Time Samples and Array Elements
ntstart = tstart*fs; % Start Time Offset in Samples
t = (ntstart:ntstart+nt0-1)/fs; z = vsound*t/2; % Time Vector & Axial Range
x = (-(nx0-1)/2:(nx0-1)/2)*pitch; % Array Element Locations

% Time-Gain Compensation to Match Delay-and-Sum Method
tgc = sqrt(t' * ones(size(x)));
RF = RF.*tgc; % tgc: time-gain compensation factor to match delay-and-sum

% Zero-padding: Picking Next Power of 2 to Optimize Sampling and Speed
nt = 2^(nextpow2(nt0)+1+upsamp(2)); % Long Enough to Avoid Range Aliasing
nx = 2^(nextpow2(nx0)+2+upsamp(1)); % Conv Tx and Rx Apertures: Lateral Span is 2*nx0 

% Exploding Reflector Model Velocity
c = vsound; % propagation velocity (m/s)
ERMv = c/2; % adaptation to monostatic synthetic aperture

% FFT
fftRF = fftshift(fft2(RF, nt, nx)); % Go Into K-Space
f = (-nt/2:nt/2-1)*fs/nt; % Get Temporal Frequency 
kx = (-nx/2:nx/2-1)/pitch/nx; % Get Lateral Spatial Frequency 
[Kx, F] = meshgrid(kx, f); % Create f-kx mesh

% Compensate for Depth Start
fftRF = fftRF.*exp(-1i*2*pi*tstart*F);

% Linear Interpolation 
FKz = ERMv*sign(F).*sqrt(Kx.^2+F.^2/ERMv^2);
fftRF = interp2(Kx,F,fftRF,Kx,FKz,'linear',0);

% Jacobian
kz = (-nt/2:nt/2-1)'/ERMv/fs/nt;
fftRF = bsxfun(@times,fftRF,kz)./(FKz+eps);

% Uncompensate for Depth Start
fftRF = fftRF.*exp(1i*2*pi*tstart*F);

% IFFT & Migrated RF
migRF = ifft2(ifftshift(fftRF), 'symmetric');
migRF = migRF(1:nt0, 1:nx0); 

end

