clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));
addpath(genpath([pwd, '/../Functions_BeamformFastFourier']));

% Load File and Set Imaging Grid
load FieldII_PW.mat; % Cyst and Lesions Phantom

% Imaging Parameters
dBrange = [-80, 0];
c = 1540; % Sound Speed [m/s]
param.fs = 1/mean(diff(time));
param.pitch = mean(diff(rxAptPos(:,1)));
param.t0 = time(1);

% Plane-Wave Reconstruction at this Angle
param.TXangle = pw_ang_deg*pi/180;
tic; [migSIG, paramOut] = lumig(rcvdata, param); toc;
img_h = hilbert(migSIG);

% Show Images
imagesc(1000*paramOut.x, 1000*paramOut.z, ...
    db(abs(img_h))-db(max(abs(img_h(:)))), dBrange);
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
title('PW Image After Compounding'); 
axis image; colormap(gray); colorbar(); drawnow;
ylim([4, 36]); % y limits [mm]