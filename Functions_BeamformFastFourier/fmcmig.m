function [migSIG,param] = fmcmig(SIG,param)

%FMCMIG   full-matrix capture (FMC) Stolt's migration
%   MIGSIG = FMCMIG(SIG,PARAM) performs a f-k migration of the signals
%   stored in the array SIG. MIGSIG contains the migrated signals. PARAM is
%   a structure that contains all the parameter values required for the
%   migration (see below for more details).
%
%   [MIGSIG,PARAM] = FMCMIG(SIG,PARAM) also returns the complete list of
%   parameters:
%   PARAM.x and PARAM.z correspond to the coordinates at which the migrated
%   signals are returned. The x-axis is parallel to the transducer and
%   pointing from element #1 to element #N (PARAM.x = 0 at the CENTER of
%   the transducer). The z-axis is PERPENDICULAR to the transducer and
%   pointing downward (PARAM.z = 0 at the level of the transducer).
%
%   Important details on FKMIG:
%   --------------------------
%   1) The signals - typically RF signals - in SIG must be acquired using a
%      FULL MATRIX CAPTURE configuration with a linear array as used in 
%      multistatic synthetic aperture ultrasound imaging. 
%   2) A bandpass filter can be used to smooth the migrated signals (see
%      PARAM.passband below).
%
%   PARAM is a structure that contains the following fields:
%   -------------------------------------------------------
%   1) PARAM.fs: sample frequency (in Hz, REQUIRED)
%   2) PARAM.pitch: pitch of the linear transducer (in m, REQUIRED)
%   3) PARAM.c: longitudinal velocity (in m/s, default = 1540 m/s)
%   4) PARAM.t0: acquisition start time (in s, default = 0)
%   5) PARAM.virtSrcDepth: depth of virtual sources (in m, default = 0).
%            Assumes all virtual sources located at the same depth.
%            Lateral spacing between virtual sources is the same as the 
%            spacing between receive elements.
%   6) PARAM.passband: Passband of the bandpass filter (default = [0,1]).
%            Must be a 2-element vector [w1,w2], where 0<w1<w2<1.0, with
%            1.0 corresponding to half the sample rate (i.e. param.fs/2).
%            A 5th order bandpass Butterworth filter is used.
%            Example: If you measured the RF signals at a sample rate of 20
%            MHz with a linear array whose passband is 3-7 MHz, you may use
%            PARAM.passband = [3e6,7e6]/(20e6/2) = [0.3 0.7].
%
%   References
%   ---------- 
%   Cheng and Lu, Extended High-Frame Rate Imaging Method with 
%   Limited-Diffraction Beams. 
%
%   Hunter et al., The Wavenumber Algorithm for Full-Matrix Imaging Using
%   an Ultrasonic Array. 
%
%   Jiang et al., The Spectrum-Beamformer for Conventional B-Mode
%   Ultrasound Imaging System: Principle, Validation, and Robustness. 
%   
%   -- Damien Garcia & Louis Le Tarnec -- 2011/08, revised 2013/11/18
%   -- Rehman Ali -- updated 2023/02

[nt,nx,nxx] = size(SIG);
assert(nx==nxx, 'last two dimensions of SIG must be the same');

%----- Input parameters -----
%-- 1) Speed of sound
if ~isfield(param,'c')
    param.c = 1540; % longitudinal velocity in m/s
end
c = param.c;
%-- 2) Sample frequency
if ~isfield(param,'fs')
    error(['A sample frequency (fs) ',...
        'must be specified as a structure field.'])
end
%-- 3) Acquisition start time
if ~isfield(param,'t0') % in s
    param.t0 = 0; % acquisition time start in s
else
    assert(isscalar(param.t0),'PARAM.t0 must be a scalar.')
end
%-- 4) Pitch
if ~isfield(param,'pitch') % in m
    error('A pitch value (PARAM.pitch) is required.')
end
%-- 5) Virtual Source Depth
if ~isfield(param,'virtSrcDepth')
    param.virtSrcDepth = 0;
end
assert(isscalar(param.virtSrcDepth),...
    'PARAM.virtSrcDepth must be a single scalar value')
%-- 6) Passband filter
if isfield(param,'passband')
    wn = param.passband;
else
    wn = [0 1];
end
assert(numel(wn)==2 & wn(1)<wn(2) & wn(1)>=0 & wn(2)<=1 ,...
    'PARAM.passband must be a 2-element vector [w1,w2] with 0<w1<w2<1.')
%----- end of Input parameters -----

SIG = double(SIG);

%-- Apply a bandpass filter if a passband is given
if ~isequal(wn(:),[0;1])
    [b,a] = butter(5,wn);
    SIG = filtfilt(b,a,SIG);
end

%-- Temporal shift
t0 = param.t0;
ntshift = max(round(t0*param.fs));

%-- Zero-padding before FFTs
%- time direction
ntFFT = 2*nt+ntshift;
if rem(ntFFT,2)==1 % ntFFT must be even
    ntFFT = ntFFT+1;
end
%- x-direction
factor = 1.5;
nxFFT = ceil(factor*nx); % in order to avoid lateral edge effects
if rem(nxFFT,2)==1 % nxFFT must be even
    nxFFT = nxFFT+1;
end

%-- Frequency axes
f0 = (-ntFFT/2+1:ntFFT/2)'*param.fs/ntFFT; 
k0 = f0/c; kz = 2*k0;
kuv = (-nxFFT/2+1:nxFFT/2)/param.pitch/nxFFT;
[ku,k,kv] = meshgrid(kuv,k0,kuv); 
validity_region = (k.^2>ku.^2) & (k.^2>kv.^2); %- Validity Regions
krz = sqrt(k.^2-ku.^2); ktz = sqrt(k.^2-kv.^2); 
directivity = krz .* ktz; %- Directivity Weighting
virtSrcShift = ... %- Virtual Source Shifting Factors
    exp(1i*2*pi*(k-ktz)*param.virtSrcDepth); 
startTimeShift = exp(-1i*2*pi*k*c*t0); %- Start Time Shift

%-- FFT each axis
SIG = fftshift(fftn(SIG,[ntFFT,nxFFT,nxFFT]));
%-- Apply directivity weightings and start time corrections
SIG = directivity .* virtSrcShift .* startTimeShift .* SIG;

%-- FMC Stolt Migration
migSIG = zeros(ntFFT,2*nxFFT-1);
[Kv,Kz] = meshgrid(kuv,kz);
for ku_idx = 1:nxFFT
    % Extract Single Slice in Ku
    SIGku = squeeze(SIG(:,ku_idx,:));
    % Stolt Mapping
    K = real(sqrt(Kz.^4+2*(kuv(ku_idx)^2+Kv.^2).*(Kz.^2)+...
        (kuv(ku_idx)^4)+(Kv.^4)-2*(kuv(ku_idx)^2).*(Kv.^2)))./(2*Kz);
    valid_k = all(squeeze(validity_region(:,ku_idx,:)),2);
    migSIGku = zeros(ntFFT,nxFFT);
    for kv_idx = 1:nxFFT
        migSIGku(:,kv_idx) = interp1(k0(k0(:)>0 & valid_k), ...
            SIGku(k0(:)>0 & valid_k, kv_idx), K(:,kv_idx), 'spline', 0);
    end
    % Exclusion Zone: Remove Incorrectly Mapped Signals
    migSIGku(((Kv.^2+Kz.^2)<=((kuv(ku_idx)).^2)) | ...
        ((Kv.^2-Kz.^2)>=(kuv(ku_idx)).^2) | (Kz == 0)) = 0;
    % Stolt Map this Slice and Accumulate
    migSIG(:, (1:nxFFT)+(ku_idx-1)) = ...
        migSIG(:, (1:nxFFT)+(ku_idx-1)) + migSIGku;
end

%-- f-IFFT
[~, f] = meshgrid(1:2*nxFFT-1,f0);
depth_start = exp(1i*2*pi*t0*f); %- Modulation for depth start
migSIG = ifftn(ifftshift(migSIG.*depth_start)); %- IFFT
migSIG = migSIG(1:nt, 1:2*nx-1); %-- Final migrated signal

%-- Grid coordinates
if nargin>1
    param.x = (-(nx-1):(nx-1))*(param.pitch/2);
    ntstart = t0*param.fs; % Start Time Sample
    t = (ntstart:ntstart+nt-1)/param.fs; 
    param.z = c*t/2; % Time Vector & Axial Range
end

end