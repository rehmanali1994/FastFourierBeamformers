clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));
addpath(genpath([pwd, '/../Functions_BeamformDAS']));

% Load File Saved by GenFullSynthData.m
load('../Examples_FMC/FieldII.mat'); % Cyst and Lesions Phantom
scat = scat/max(abs(scat(:)));

% Apply Bandpass Filter to Channel Data
passbandHz = [1e6,10e6]; % Passband [Hz]
wn = passbandHz/(fs/2);
[b,a] = butter(5,wn);
scat = filtfilt(b,a,scat);

% Setup Transmit Imaging Case Here
txAptPos = rxAptPos; % Using Same Array to Transmit and Recieve
c = 1540; % Speed of Sound [m/s] for Transmit Beamforming
pw_ang_deg = linspace(-30,30,61); % Plane-Wave Angle [Degrees]

% Get Different Rx Data for Each Tx 
rcvdata = zeros(numel(time), no_elements, numel(pw_ang_deg));
parfor ang_idx = 1:numel(pw_ang_deg)
    % Transmit Direction
    tx_dir = [-sind(pw_ang_deg(ang_idx)),0,cosd(pw_ang_deg(ang_idx))]; 
    % Transmit Center
    if pw_ang_deg(ang_idx) > 0
        tx_center = txAptPos(1,:);
    else
        tx_center = txAptPos(end,:);
    end
    % Convert Full-Matrix Capture (FMC) Dataset to Plane-Wave (PW) Dataset
    rcvdata(:,:,ang_idx) = ...
        focus_fs_to_TxBeam(time, scat, ...
        rxAptPos, txAptPos, tx_center, ...
        tx_dir, Inf, ones(1,no_elements), 0, c);
    disp(['Completed Transmit ', num2str(ang_idx), ...
        '/', num2str(numel(pw_ang_deg))]);
end

% Save Focused Transmit Data to File
save('FieldII_PW.mat','-v7.3','time',...
    'rxAptPos','txAptPos','pw_ang_deg','rcvdata');