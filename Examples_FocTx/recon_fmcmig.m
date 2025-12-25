clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));
addpath(genpath([pwd, '/../Functions_BeamformFastFourier']));

% Load File and Set Imaging Grid
load FieldII_FocTx.mat; % Cyst and Lesions Phantom

% Imaging Parameters
dBrange = [-80, 0];
c = 1540; % Sound Speed [m/s]
param.fs = 1/mean(diff(time));
param.pitch = mean(diff(rxAptPos(:,1)));
param.t0 = time(1);
param.virtSrcDepth = tx_focDepth;

% Plane-Wave Reconstruction at this Angle
tic; [migSIG, paramOut] = fmcmig(rcvdata, param); toc;

% Show Images
imagesc(1000*paramOut.x, 1000*paramOut.z, ...
    db(abs(migSIG))-db(max(abs(migSIG(:)))), dBrange);
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
title('Focused TX Image After Compounding'); 
axis image; colormap(gray); colorbar(); drawnow;
ylim([4, 36]); % y limits [mm]