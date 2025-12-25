clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));
addpath(genpath([pwd, '/../Functions_BeamformDAS']));

% Load File and Set Imaging Grid
load FieldII_PW.mat; % Cyst and Lesions Phantom

% Imaging Parameters
dBrange = [-60,0];
c = 1540; % Sound Speed [m/s]
param.fs = 1/mean(diff(time));
param.pitch = mean(diff(rxAptPos(:,1)));
param.t0 = time(1);
param.fnumber = 1; % Receive F-Number

% Reconstruct Image
img_h_cmpd = zeros(numel(time),size(rxAptPos,1));
for ang_idx = 1:numel(pw_ang_deg)
    % Plane-Wave Reconstruction at this Angle
    param.TXangle = pw_ang_deg(ang_idx)*pi/180;
    [migSIG, paramOut] = pw_dasmig(rcvdata(:,:,ang_idx), param);
    img_h = hilbert(migSIG);
    % Mask to Remove Anything Outside the Support of the Plane Wave
    [X_IMG, Z_IMG] = meshgrid(paramOut.x, paramOut.z);
    MASK = (X_IMG<rxAptPos(end,1)+Z_IMG*tand(pw_ang_deg(ang_idx))) & ...
        (X_IMG>rxAptPos(1,1)+Z_IMG*tand(pw_ang_deg(ang_idx)));
    % Coherently Compound Image
    img_h_cmpd = img_h_cmpd + MASK .* img_h;
    % Show Images
    subplot(1,4,1); imagesc(1000*paramOut.x, 1000*paramOut.z, ...
        db(abs(img_h))-db(max(abs(img_h(:)))), dBrange);
    xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
    title('PW Image Before Masking');
    axis image; colormap(gray); colorbar();
    subplot(1,4,2); imagesc(1000*paramOut.x, 1000*paramOut.z, MASK);
    xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
    title('PW Mask'); axis image; colormap(gray); colorbar();
    subplot(1,4,3); imagesc(1000*paramOut.x, 1000*paramOut.z, ...
        db(abs(MASK.*img_h))-db(max(abs(img_h(:)))), dBrange);
    xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
    title('PW Image After Masking'); 
    axis image; colormap(gray); colorbar(); 
    subplot(1,4,4); imagesc(1000*paramOut.x, 1000*paramOut.z, ...
        db(abs(img_h_cmpd))-db(max(abs(img_h_cmpd(:)))), dBrange);
    xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
    title('PW Image After Compounding'); 
    axis image; colormap(gray); colorbar(); drawnow;
end
