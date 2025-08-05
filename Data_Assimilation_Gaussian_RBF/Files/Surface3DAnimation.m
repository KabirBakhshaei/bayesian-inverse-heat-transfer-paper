% 3D Surface Plot Animation
% This MATLAB script generates a 3D surface plot animation to visualize
% the true and reconstructed heat flux over time. It loads data from
% specific files, calculates the heat flux reconstruction, and creates an
% animation showing the true and reconstructed heat flux in a 3D plot.
% The animation is saved as '3D Combined surface plot.avi' and '3D Combined surface plot.gif'.
%
% Please make sure to set the appropriate file paths for data loading.
clc; clearvars; close all;

outputDir = '../Results';  % or full path
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
aviFilename = fullfile(outputDir, '3D Combined surface plot.avi');
gifFilename = fullfile(outputDir, '3D Combined surface plot.gif');

% Load the data for the first animation
gtrue = load('./ITHACAoutput/projection/TrueHeatFlux/HeatFluxTrue_mat.txt'); % [101, 400]
gtrue = gtrue(1:end-1,:);

% Load parameterMean and heatFluxSpaceRBF matrices
parameterMean = load('./ITHACAoutput/reconstruction/parameterMean_mat.txt'); % [25, 100]
heatFluxSpaceRBF = load('./ITHACAoutput/projection/HeatFluxSpaceRBF/heat_flux_space_basis_mat.txt'); % [25, 400]

% Get the dimensions of parameterMean and heatFluxSpaceRBF matrices
[n1, m1] = size(parameterMean); % n1 = 25 (mean/weight), m1 = 100 (times)
[n, m] = size(heatFluxSpaceRBF); % n = 25 (RBF), m = 400 (faces)

% Initialize the out matrix with zeros
out = zeros(m1, m); % (100, 400)

% Perform matrix multiplication and summation
for i = 1:m1 % Loop over times (100 times)
    for j = 1:m % Loop over faces (400 faces)
        out(i, j) = sum(parameterMean(:, i) .* heatFluxSpaceRBF(:, j)); % Matrix multiplication and summation
    end
end

% Load the 'xyz.txt' file or prompt user to run the notebook if not found
xyz_txt_path = fullfile(outputDir, 'xyz.txt');

if ~isfile(xyz_txt_path)
    error(['The file "xyz.txt" was not found in ', outputDir, ...
           '. Please run "plots.ipynb" first to generate it from "xyz.npy".']);
end
% Load the 'xyz.npy' file
xyz = load(xyz_txt_path);

% Assume x, y, z are reshaped properly or use griddata if needed
x = xyz(1:3:end); % Assume these are grid coordinates or use meshgrid
y = xyz(2:3:end);
z = xyz(3:3:end);

% Load trueTimeVec_mat.txt
timeInstants = load('./ITHACAoutput/true/trueTimeVec_mat.txt');

% Create a VideoWriter object to save the animation in AVI format
aviObj = VideoWriter(aviFilename, 'Motion JPEG AVI');
aviObj.FrameRate = 10;  % Adjust the frame rate as needed
aviObj.Quality = 100;   % Set the video quality
open(aviObj);

% Create a figure for the 3D plots and set it to full screen
fig = figure('Position', get(0, 'ScreenSize'));

xmin = 0;
xmax = 2;
ymin = 0;
ymax = 1;
zmin = 0;
zmax = 1.2;
gtrue_min = min(min(gtrue(:)), min(out(:)));
gtrue_max = max(max(gtrue(:)), max(out(:)));

% Define the grid for surface plot
[X, Y] = meshgrid(linspace(xmin, xmax, sqrt(m)), linspace(ymin, ymax, sqrt(m)));

% Initialize variables for GIF
delayTime = 0.2; % Delay time in seconds

for i = 1:length(timeInstants)
    % Reshape gtrue and out for surface plot
    gtrue_t = reshape(gtrue(i, :), sqrt(m), sqrt(m));
    out_t = reshape(out(i, :), sqrt(m), sqrt(m));

    % Plot true heat flux
    subplot(1, 2, 1);
    surf(X, Y, gtrue_t, 'EdgeColor', 'none');
    colormap('jet');
    xlabel('X(m)', 'FontWeight', 'bold');
    ylabel('Z(m)', 'FontWeight', 'bold');
    zlabel('True Heat Flux (W/m^2)', 'FontWeight', 'bold');
    %title(['True Heat Flux at Time Instant ', num2str(timeInstants(i))]);
    title(['True Heat Flux at Time Instant ', num2str(timeInstants(i), '%.1f'), ' seconds']);
    axis([xmin xmax ymin ymax gtrue_min gtrue_max]);
    view(3); % Sets the view to 3D perspective

    % Plot reconstructed heat flux
    subplot(1, 2, 2);
    surf(X, Y, out_t, 'EdgeColor', 'none');
    colormap('jet')
    xlabel('X(m)', 'FontWeight', 'bold');
    ylabel('Z(m)', 'FontWeight', 'bold');
    zlabel('Reconstructed Heat Flux (W/m^2)', 'FontWeight', 'bold');
    %title(['Reconstructed Heat Flux at Time Instant ', num2str(i)]);
    title(['Reconstructed Heat Flux at Time Instant ', num2str(timeInstants(i), '%.1f'), ' seconds using Gaussian RBF']);


    axis([xmin xmax ymin ymax gtrue_min gtrue_max]);
    view(3);

    % Capture the frame and add it to the AVI file
    frame = getframe(fig);
    writeVideo(aviObj, frame);

    % Add frame to GIF file
    [imind, cm] = rgb2ind(frame.cdata, 64);
    if i == 1
        imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', delayTime);
    else
        imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end

    % Optional pause
    %pause(0.1);
end

% Close the video file
close(aviObj);

% Display a message indicating that the animations have been saved
disp(['GIF will be saved to: ', gifFilename]);
disp(['AVI will be saved to: ', aviFilename]);

