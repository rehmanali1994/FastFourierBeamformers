clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));
addpath(genpath([pwd, '/../Functions_BeamformDAS']));

% Load File and Set Imaging Grid
load FieldII.mat; % Cyst and Lesions Phantom
num_x = 255; xlims = (12.7e-3)*[-1, 1];
num_z = 1326; zlims = [4e-3, 36e-3];

% Select Subset of Transmit Elements
tx_elmts = 1:128;
txAptPos = rxAptPos(tx_elmts,:);
rxdata_h = hilbert(scat(:,:,tx_elmts));

% Points to Focus and Get Image At
x_img = linspace(xlims(1), xlims(2), num_x);
z_img = linspace(zlims(1), zlims(2), num_z);
dBrange = [-80, 0]; c = 1540;

% Full Synthetic Aperture Image Reconstruction
[X, Y, Z] = meshgrid(x_img, 0, z_img);
foc_pts = [X(:), Y(:), Z(:)];
tic; focData = bfm_fs_fast(time, rxdata_h, foc_pts, rxAptPos, txAptPos, 0, 0, c); toc;
img_h = reshape(focData, [numel(x_img), numel(z_img)])';
imagesc(1000*x_img, 1000*z_img, 20*log10(abs(img_h)/max(abs(img_h(:)))), dBrange); 
axis image; xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title('DAS Beamforming'); colormap(gray); colorbar();