% -------------------------------------------------------------------------
% MATLAB Code for Creating PNG images of Recunstrcuted and True Heat Flux
% for the paper
% -------------------------------------------------------------------------
%
% This MATLAB code demonstrates how to create PNG images of a series of
% heat flux contour plots.
%
% Instructions:
% 1. Make sure you have the necessary data files in the specified paths.
% 2. Adjust video settings, such as frame rate and file name, as needed.
% 3. Run the code to generate the AVI video.
% -------------------------------------------------------------------------

clc; clearvars; close all;

%% Laod the true HF
gtrue = load('./ITHACAoutput/projection/TrueHeatFlux/HeatFluxTrue_mat.txt'); % [101, 400]
gtrue = gtrue (1:end-1,:);

%% Creating reconstructed HF
% Load parameterMean and heatFluxSpaceRBF matrices
parameterMean = load('./ITHACAoutput/reconstruction/parameterMean_mat.txt'); % [5, 100]
heatFluxSpaceRBF = load('./ITHACAoutput/projection/HeatFluxSpaceRBF/heat_flux_space_basis_mat.txt'); % [5, 400]

% Get the dimensions of parameterMean and heatFluxSpaceRBF matrices
[n1, m1] = size(parameterMean);    % n1 = 5 (mean/weight), m1 = 100 (times)
[n, m] = size(heatFluxSpaceRBF);   % n = 5 (RBF),          m = 400 (faces)

% Initialize the out matrix with zeros
out = zeros(m1, m); % (100, 400)

% Perform matrix multiplication and summation
for i = 1:m1       % Loop over times (100 times)
    for j = 1:m    % Loop over faces (400 faces)
        out(i, j) = sum(parameterMean(:, i) .* heatFluxSpaceRBF(:, j)); % Matrix multiplication and summation
    end
end
%% Load time vector
% Load trueTimeVec_mat.txt
timeInstants = load('./ITHACAoutput/true/trueTimeVec_mat.txt');

%% Load the 'xyz.npy' file
outputDir = '../Results';  % or full path

% Load the 'xyz.txt' file or prompt user to run the notebook if not found
xyz_txt_path = fullfile(outputDir, 'xyz.txt');

if ~isfile(xyz_txt_path)
    error(['The file "xyz.txt" was not found in ', outputDir, ...
           '. Please run "plots.ipynb" first to generate it from "xyz.npy".']);
end
% Load the 'xyz.npy' file
xyz = load(xyz_txt_path);


% Extract x, y, and z arrays
X_coordinate = xyz(1:3:end);   % Start at index 1 and take every third element
y_coordinate = zeros(m,1);     %xyz(2:3:end);  % Start at index 2 and take every third element
z_coordinate = xyz(3:3:end);   % Start at index 3 and take every third element

%% Meshgrid       Create a grid of points that matches the dimensions of heatFlux
%[Xgrid, Zgrid] = meshgrid(X_coordinate,z_coordinate);
fineGridResolution = 0.1; % Define a finer grid resolution
%[Xgrid, Zgrid] = meshgrid(min(X_coordinate):fineGridResolution:max(X_coordinate), min(z_coordinate):fineGridResolution:max(z_coordinate));
[Xgrid, Zgrid] = meshgrid(0:fineGridResolution:2, 0:fineGridResolution:1.2);
Ygrid = zeros(size(Xgrid)); % Set Y-coordinate to zero for all points

%% 
%gTrue_gRecunstructed_min = min(min(gtrue(:)), min(out(:)));
%gTrue_gRecunstructed_max = max(max(gtrue(:)), max(out(:)));

gTrue_gRecunstructed_min = min(gtrue(:));
gTrue_gRecunstructed_max = max(gtrue(:));

fontSize = 12;                % Set the desired font size
fontType = 'Times New Roman'; % Set the desired font type
fontSize2 = 10;

for i = 1: m1

    figure(1)
    subplot(2,1,1);
    heatFlux = out(i,:);

    % Interpolate the reconstructed heat flux values to match the grid
    heatFluxInterpolated = griddata(X_coordinate, z_coordinate, heatFlux, Xgrid, Zgrid, 'nearest'); %cubic, nearest:  Use cubic interpolation for smoother contours

    % Compute the minimum and maximum values of reconstructed heat flux for this entire time steps
    minHeatFluxReconstructed = gTrue_gRecunstructed_min; %min(out(:)); % min(heatFluxInterpolated(:));
    maxHeatFluxReconstructed = gTrue_gRecunstructed_max; %max(out(:)); % max(heatFluxInterpolated(:)); 

    % Specifying the number of contour levels and their values
    numLevels = 100;  % Creating a contour plot with 10 contour lines representing different heat flux levels.
    contourLevels = linspace(minHeatFluxReconstructed, maxHeatFluxReconstructed, numLevels);

    % Create the contour plot
    %contourHandle = contourf(Xgrid, Zgrid, heatFluxInterpolated, contourLevels, 'LineWidth', 5);
    contourHandleReconstructed = contourf(Xgrid, Zgrid, heatFluxInterpolated, contourLevels, 'LineStyle', 'none');

    % Adding labels, titles, and color map to make the contour plot more informative:
    %xlabel('X');
    %ylabel('Z');

    xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
    ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');

    ax = gca;
    ax.XAxis.FontSize = fontSize2; ax.XAxis.FontName = fontType; ax.XAxis.FontWeight = 'bold'; 
    ax.YAxis.FontSize = fontSize2; ax.YAxis.FontName = fontType; ax.YAxis.FontWeight = 'bold';

    grid on
    title(['Transparent Cube with Contour Plot of Reconstructed Heat Flux at Time Instant ', num2str(i)]);
    CC = colorbar; % Add a colorbar to the plot
    CC.Label.String = 'Reconstructed HF (W/m^{2})';
    CC.Label.FontSize = 14;
    colormap('jet'); % parula
    clim([minHeatFluxReconstructed maxHeatFluxReconstructed])
    axis equal
        %%
        % Check if i is one of the specified values (25, 50, 75, 100)
    if ismember(i, [25, 50, 75, 100])
        % Create a new subplot
        figure(2); % You may want to use a different figure number
        subplot(4,1,find(i == [25, 50, 75, 100]));
        
        contourHandleReconstructed = contourf(Xgrid, Zgrid, heatFluxInterpolated, contourLevels, 'LineStyle', 'none');

        % Customize labels, titles, etc. for the additional subplot
        xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
        ylabel('Z','FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
        ax = gca;
        ax.XAxis.FontSize = fontSize2; ax.XAxis.FontName = fontType; ax.XAxis.FontWeight = 'bold'; 
        ax.YAxis.FontSize = fontSize2; ax.YAxis.FontName = fontType; ax.YAxis.FontWeight = 'bold';

        title(['Time = ', num2str(i * 0.2), ' sec'], 'FontSize', 10);
        CC = colorbar; % Add a colorbar to the plot
        %CC.Label.String = 'Reconstructed HF (W/m^{2})';
        CC.Label.FontSize = 10;
        colormap('jet'); % parula
        clim([minHeatFluxReconstructed maxHeatFluxReconstructed])
        axis equal
    end
    fig = gcf;
    set(fig, 'WindowState', 'maximized');
%%
    figure(1)
    subplot(2,1,2);
    heatFluxTrue = gtrue(i,:);

    % Interpolate the heat flux values to match the grid
    heatFluxInterpolatedTrue = griddata(X_coordinate, z_coordinate, heatFluxTrue, Xgrid, Zgrid, 'nearest'); % cubic, nearest:  Use cubic interpolation for smoother contours
 
    % Compute the minimum and maximum values of (out/gtrue/heatFluxInterpolated) for this (time instant/entire)
    minHeatFluxTrue = gTrue_gRecunstructed_min; %min(gtrue(:)); % min(heatFluxInterpolatedTrue(:));
    maxHeatFluxTrue = gTrue_gRecunstructed_max; %max(gtrue(:)); % max(heatFluxInterpolatedTrue(:)); 

    % Specifying the number of contour levels and their values
    numLevelsTrue = 100;  % Creating a contour plot with 10 contour lines representing different heat flux levels.
    contourLevelsTrue = linspace(minHeatFluxTrue, maxHeatFluxTrue, numLevelsTrue);

    % Create the contour plot
    %contourHandle = contourf(Xgrid, Zgrid, heatFluxInterpolated, contourLevels, 'LineWidth', 5);
    contourHandleTrue = contourf(Xgrid, Zgrid, heatFluxInterpolatedTrue, contourLevelsTrue, 'LineStyle', 'none');

    % Adding labels, titles, and color map to make the contour plot more informative:
    %xlabel('X');
    %ylabel('Z');

    xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
    ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');

    ax = gca;
    ax.XAxis.FontSize = fontSize2; ax.XAxis.FontName = fontType; ax.XAxis.FontWeight = 'bold'; 
    ax.YAxis.FontSize = fontSize2; ax.YAxis.FontName = fontType; ax.YAxis.FontWeight = 'bold';

    grid on
    title(['Transparent Cube with Contour Plot of True Heat Flux at Time Instant ', num2str(i)]);
    CC = colorbar; % Add a colorbar to the plot
    CC.Label.String = 'True HF (W/m^{2})';
    CC.Label.FontSize = 14;
    colormap('jet'); % parula
    clim([minHeatFluxTrue maxHeatFluxTrue])
    axis equal

    %%
    if ismember(i, [25, 50, 75, 100])
        % Create a new subplot
        figure(3); % You may want to use a different figure number
        subplot(4,1,find(i == [25, 50, 75, 100]));
        
        contourHandleTrue = contourf(Xgrid, Zgrid, heatFluxInterpolatedTrue, contourLevelsTrue, 'LineStyle', 'none');

        % Customize labels, titles, etc. for the additional subplot
        xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
        ylabel('Z','FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
        ax = gca;
        ax.XAxis.FontSize = fontSize2; ax.XAxis.FontName = fontType; ax.XAxis.FontWeight = 'bold'; 
        ax.YAxis.FontSize = fontSize2; ax.YAxis.FontName = fontType; ax.YAxis.FontWeight = 'bold';

        title(['Time = ', num2str(i * 0.2), ' sec'], 'FontSize', 10);
        CC = colorbar; % Add a colorbar to the plot
        %CC.Label.String = 'Reconstructed HF (W/m^{2})';
        CC.Label.FontSize = 10;
        colormap('jet'); % parula
        clim([minHeatFluxReconstructed maxHeatFluxReconstructed])
        axis equal

        fig = gcf;
        set(fig, 'WindowState', 'maximized');
        
                % ----------- New Figure 4: Relative Error (%) -----------
        figure(4);
        subplot(4,1,find(i == [25, 50, 75, 100]));

        % Calculate relative error in percentage
        epsilon = 1e-8;
        relativeError = abs(heatFluxTrue - heatFlux) ./ (abs(heatFluxTrue) + epsilon) * 100;
        relativeErrorInterpolated = griddata(X_coordinate, z_coordinate, relativeError, Xgrid, Zgrid, 'nearest');

        % Plot the relative error
        contourf(Xgrid, Zgrid, relativeErrorInterpolated, 100, 'LineStyle', 'none');
        xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
        ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');

        ax = gca;
        ax.XAxis.FontSize = fontSize2; ax.XAxis.FontName = fontType; ax.XAxis.FontWeight = 'bold';
        ax.YAxis.FontSize = fontSize2; ax.YAxis.FontName = fontType; ax.YAxis.FontWeight = 'bold';

        title(['Time = ', num2str(i * 0.2), ' sec'], 'FontSize', 10);
        CC = colorbar;
        CC.Label.String = 'Relative Error (%)';
        CC.Label.FontSize = 10;
        colormap('jet');
        caxis([0 100]);  % adjust range if needed
        axis equal

        fig = gcf;
        set(fig, 'WindowState', 'maximized');

    end
    %%
    % Get the handle to the current figure
    fig = gcf;

    % Set the figure's position and size to make it full screen
    %set(fig, 'Position', [0, 0, 1920, 1080]); % Adjust the values as needed
    % Maximize the figure
    set(fig, 'WindowState', 'maximized');

    % Rotate the view for a better perspective
    % view(0,90);
    %view(3)
    %view([1,0.6,0.05])

    % Pause to allow time to view the plot (optional)
    pause(0.01);

%     if i == 5
%         break
%     end

    %%

    % Set the figure's position and size to make it full screen
    %set(fig, 'Position', [0, 0, 1920, 1080]); % Adjust the values as needed
    % Maximize the figure
    set(fig, 'WindowState', 'maximized');
end

% Define filenames
filenameFigure1 = 'ComparisnCountors.png';
filenameFigure2 = 'SnapshotCountorsMultiquadric.png';
filenameFigure3 = 'SnapshotCountorsTrue.png';

% Delete existing files if they exist
if exist(filenameFigure1, 'file')
    delete(filenameFigure1);
end

if exist(filenameFigure2, 'file')
    delete(filenameFigure2);
end

if exist(filenameFigure3, 'file')
    delete(filenameFigure3);
end

%saveas(figure(1), 'ComparisnCountors.png', 'png');
%saveas(figure(2), 'SnapshotCountorsReconstructed.png', 'png');
%saveas(figure(3), 'SnapshotCountorsTrue.png', 'png', 'Resolution');

%print(figure(1), fullfile(outputDir, 'ComparisnCountors.png'), '-dpng', '-r400'); % Adjust '-r' value as needed for resolution
print(figure(2), fullfile(outputDir, 'Figure 14b SnapshotCountorsMultiquadric.png'), '-dpng', '-r400');
print(figure(3), fullfile(outputDir, 'Figure 14C SnapshotCountorsTrue.png'), '-dpng', '-r400');




filenameFigure4 = 'Figure 14e SnapshotCountorsMultiquadricRelative.png';
if exist(filenameFigure4, 'file')
    delete(filenameFigure4);
end
print(figure(4), fullfile(outputDir, filenameFigure4), '-dpng', '-r400');
