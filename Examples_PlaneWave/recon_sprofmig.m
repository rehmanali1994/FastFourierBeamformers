clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));

% Load File and Set Imaging Grid
load FieldII_PW.mat; % Cyst and Lesions Phantom

% Imaging Parameters
dBrange = [-80, 0];
c = 1540; % Sound Speed [m/s]
param.fs = 1/mean(diff(time));
param.pitch = mean(diff(rxAptPos(:,1)));
param.t0 = time(1);

% Delay and Apodization
nx = size(txAptPos,1);
param.apod = ones(numel(pw_ang_deg),nx);
param.delays = zeros(size(param.apod));
for ang_idx = 1:numel(pw_ang_deg)
    sinA = sind(pw_ang_deg(ang_idx)); % sine of steering angle
    param.delays(ang_idx,:) = ...
        -sinA*((nx-1)*(sinA<0)-(0:nx-1))*param.pitch/c; 
end

% Selecting a Subset of Transmits
tx_evts = 1:1:61;
rcvdata = rcvdata(:,:,tx_evts);
param.apod = param.apod(tx_evts,:);
param.delays = param.delays(tx_evts,:);

% Plane-Wave Reconstruction at this Angle
tic; [migSIG, paramOut] = sprofmig(rcvdata, param); toc;
% Changed inside sprofmig: factor = 2.5;
% sprofmig needs more lateral zeropadding than rfmcmig

% Show Images
imagesc(1000*paramOut.x, 1000*paramOut.z, ...
    db(abs(migSIG))-db(max(abs(migSIG(:)))), dBrange);
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
title('PW Image After Compounding'); 
axis image; colormap(gray); colorbar(); drawnow;
ylim([4, 36]); % y limits [mm]