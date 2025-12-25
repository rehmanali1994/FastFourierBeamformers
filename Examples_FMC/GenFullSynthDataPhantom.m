clear
clc

% Changing to Field II folder and running field_init
field_init(-1);

% Global Constants
c = 1540; % Speed of Sound [m/s]
fs = 25e6; % Sampling Frequency [Hz]
set_field('c',c); % [m/s]
set_field('fs',fs); % [Hz]

% Settings
fTx = 6e6; % Transmit Frequency [Hz]
lambda = c/fTx; % Wavelength [m]

% Transducer information
no_elements = 128;  % number of physical elements.
pitch = 0.2e-3; % [m] pitch or element spacing
width = pitch; % [m] width of elements in x-direction
height = 6e-3;  % [m] width of elements in y-direction
kerf = pitch-width; % [m] distance between elements (in x-direction)
no_sub_x = 1;  % number of sub-divisions of elements in x-direction
no_sub_y = 1;  % number of sub-divisions of elements in y-direction
x = (-(no_elements-1)/2:(no_elements-1)/2)*pitch; % Aperture Element Positions
rxAptPos = [x', zeros(size(x')), zeros(size(x'))]; % Defining Rx Aperture Positions

% Transmit Focusing
focDepth = 2.0e-2; % Focal Depth [m]
focus = [0,0,focDepth];  % [m] Fixed focus for array (x,y,z). Vector with three elements.

% Point Target Phantom
% Step 1 Create Some Speckle
numPoints = 40000;
point_pos_z = (4e-3) + (3.2e-2)*rand(numPoints,1);
point_pos_x = min(x)+(max(x)-min(x))*rand(numPoints,1);
point_positions = [point_pos_x, zeros(numPoints,1), point_pos_z];
amplitudes = randn(numPoints,1);
% Step 2: Add Anechoic Lesions
lesion_ctr_pos_x = (-7.5e-3):(5.0e-3):(7.5e-3);
lesion_ctr_pos_z = (7.5e-3):(5.0e-3):(32.5e-3);
for lesion_ctr_x = lesion_ctr_pos_x
    for lesion_ctr_z = lesion_ctr_pos_z
        amplitudes(sqrt((point_pos_x-lesion_ctr_x).^2 + ...
            (point_pos_z-lesion_ctr_z).^2) < (1.5e-3)) = ...
            0*amplitudes(sqrt((point_pos_x-lesion_ctr_x).^2 + ...
            (point_pos_z-lesion_ctr_z).^2) < (1.5e-3));
    end
end
% Step 3 Point Targets
numRows = 7; numCols = 5;
point_pos_z = (5e-3) + (30e-3)*(0:numRows-1)/(numRows-1);
point_pos_x = (-10e-3) + (20e-3)*(0:numCols-1)/(numCols-1);
[point_pos_X, point_pos_Z] = meshgrid(point_pos_x, point_pos_z);
point_positions = [point_positions; ...
    point_pos_X(:), zeros(numRows*numCols,1), point_pos_Z(:)];
amplitudes = [amplitudes; 30*ones(numCols*numRows,1)];

% Set Up Transducer Geometry:
Tx = xdc_linear_array(no_elements, width, height, kerf, no_sub_x, no_sub_y, focus);
Rx = xdc_linear_array(no_elements, width, height, kerf, no_sub_x, no_sub_y, focus);

% Transmit Pulse and Impulse Response
fracBW = 0.7;  % Fractional Bandwidth of the Transducer
tc = gauspuls('cutoff', fTx, fracBW, -6, -80);  % Cutoff time at -40dB, fracBW @ -6dB
t = -tc:1/fs:tc;  % (s) Time Vector centered about t=0
impResp = gauspuls(t,fTx,fracBW); % Calculate Transmit Pulse
xdc_impulse(Tx,impResp); % Transmit Gaussian Pulse
xdc_impulse(Rx,1); % (impulse)

% Gets Lateral Transmit Beam Plot at Focal Depth
[scat, start_t] = calc_scat_all(Tx, Rx, point_positions, amplitudes, 1);
time = (0:size(scat,1)-1).*(1/fs) + (start_t-tc);
scat = reshape(scat, [size(scat,1), no_elements, no_elements]);

% Save Multistatic Data
save -v7.3 FieldII.mat;
disp('Data has been saved to file');