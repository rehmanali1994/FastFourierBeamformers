clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));
addpath(genpath([pwd, '/../Functions_BeamformFastFourier']));

% Load File and Set Imaging Grid
load FieldII.mat; % Cyst and Lesions Phantom

% Select Subset of Transmit Elements
rxdata = scat; % Multistatic Dataset

% Construct Monostatic Dataset
mono_dataset = zeros(numel(time), no_elements);
for elmt_idx = 1:no_elements
    mono_dataset(:,elmt_idx) = rxdata(:,elmt_idx,elmt_idx);
end

% Imaging Parameters
dBrange = [-60, 0]; c = 1540;
RF = mono_dataset;
fs = 1./mean(diff(time));
tstart = time(1);
pitch = mean(diff(rxAptPos(:,1)));
vsound = c; % Sound Speed [m/s]
upsamp = [1,1];

% Monostatic Stolt Reconstruction
tic; [x, z, migRF] = fkmig_monostat(RF, fs, tstart, pitch, vsound, upsamp); toc;

% Show Images
imagesc(1000*x, 1000*z, ...
    db(abs(hilbert(migRF)))-max(db(abs(hilbert(migRF(:))))), dBrange);
axis image; xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title('Stolt f-k Migration'); colormap(gray); colorbar();
