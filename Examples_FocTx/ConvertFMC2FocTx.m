clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));
addpath(genpath([pwd, '/../Functions_BeamformDAS']));

% Load File Saved by GenFullSynthData.m
load('../Examples_FMC/FieldII.mat'); % Cyst and Lesions Phantom

% Setup Transmit Imaging Case Here
txAptPos = rxAptPos; % Using Same Array to Transmit and Recieve
tx_focDepth = 0.020; % Transmit Focus in [m]
tx_dir = [0,0,1]; % Each Transmit Beam is Straight Ahead
tx_apert_elmts = 33; % Size of Transmit Aperture [elements]
c = 1540; % Speed of Sound [m/s] for Transmit Beamforming

% Construct Transmit Apodization
no_elements = size(txAptPos,1); % Number of Elements
apod = zeros(no_elements, no_elements); % Apodization
element_indices = 1:no_elements; % Element Position [elements]
tx_apert_about_ctr = (-(tx_apert_elmts-1)/2):((tx_apert_elmts-1)/2);
for beam_orig_idx = element_indices
    % Apodization on a Transmit-by-Transmit Basis
    elmt_tx = beam_orig_idx + tx_apert_about_ctr;
    elmt_tx = elmt_tx((elmt_tx >= 1) & (elmt_tx <= no_elements));
    apod(beam_orig_idx,elmt_tx) = 1; 
end

% Get Different Rx Data for Each Tx 
rcvdata = zeros(numel(time), no_elements, no_elements);
txBeamOrigins = txAptPos; % Beam Origin [m]
parfor beam_orig_idx = element_indices
    rcvdata(:,:,beam_orig_idx) = ...
        focus_fs_to_TxBeam(time, scat, ...
        rxAptPos, txAptPos, txBeamOrigins(beam_orig_idx,:), ...
        tx_dir, tx_focDepth, apod(beam_orig_idx,:), 0, c);
    disp(['Completed Transmit ', num2str(beam_orig_idx), ...
        '/', num2str(no_elements)]);
end

% Save Focused Transmit Data to File
save('FieldII_FocTx.mat','-v7.3','time','rxAptPos','txAptPos',...
    'apod','txBeamOrigins','tx_dir','tx_focDepth','rcvdata');