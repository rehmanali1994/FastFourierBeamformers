clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));
addpath(genpath([pwd, '/../Functions_BeamformFastFourier']));

% Load File and Set Imaging Grid
load FieldII_FocTx.mat; % Cyst and Lesions Phantom

% Selecting Subset of Transmit Events
tx_start = 3; tx_skip = 3; tx_end = 126;
rcvdata = rcvdata(:,:,tx_start:tx_skip:tx_end);

% Imaging Parameters
dBrange = [-80, 0];
c = 1540; % Sound Speed [m/s]
param.fs = 1/mean(diff(time));
param.pitch = mean(diff(rxAptPos(:,1)));
param.t0 = time(1);

% Delay and Apodization
param.apod = apod;
param.delays = zeros(size(apod));
txFoci = txBeamOrigins; txFoci(:,3) = tx_focDepth;
for beam_idx = 1:size(txFoci,1)
    dist2focus = sqrt(sum((txAptPos-repmat(txFoci(beam_idx,:), ...
        [size(txAptPos,1),1])).^2,2));
    param.delays(beam_idx,:) = (tx_focDepth - dist2focus)/c;
end
param.apod = param.apod(tx_start:tx_skip:tx_end,:);
param.delays = param.delays(tx_start:tx_skip:tx_end,:);

% Plane-Wave Reconstruction at this Angle
tic; [migSIG, paramOut] = rfmcmig(rcvdata, param); toc;

% Show Images
imagesc(1000*paramOut.x, 1000*paramOut.z, ...
    db(abs(migSIG))-db(max(abs(migSIG(:)))), dBrange);
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
title('Focused TX Image After Compounding'); 
axis image; colormap(gray); colorbar(); drawnow;
ylim([4, 36]); % y limits [mm]