clc; clearvars; close all;

outputDir = fullfile('..', 'Results', filesep);


% Create the directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Load the 'xyz.txt' and 'Temp.txt' file
xyz = load('./xyz.txt');

RBFsTempG0_5 = load('./TempG0_5.txt');
RBFsTempG1   = load('./TempG1.txt');
RBFsTempG2   = load('./TempG2.txt');
RBFsTempG2_5 = load('./TempG2_5.txt');

RBFsTempM0_5 = load('./TempM0_5.txt');
RBFsTempM1   = load('./TempM1.txt');
RBFsTempM3   = load('./TempM3.txt');
RBFsTempM7_5 = load('./TempM7_5.txt');

% Extract x, y, and z arrays
X_coordinate = xyz(1:3:end);   %                Start at index 1 and take every third element
m= length(X_coordinate);
Y_coordinate = zeros(m,1);     % xyz(2:3:end);  Start at index 2 and take every third element
Z_coordinate = xyz(3:3:end);   %                Start at index 3 and take every third element

%scatter(X_coordinate, Z_coordinate, markerSize=20, markerColor='k', 'filled', 'LineWidth', lineWidth=2);
scatter(X_coordinate, Z_coordinate, 22, 'k', 'filled', 'LineWidth', 2);
xlabel('X', 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Z', 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
set(gca, 'FontSize', 10);  
set(gca, 'FontWeight', 'bold');
% dpi = 400; print('XZCentersOfFaces.png', '-dpng', ['-r' num2str(dpi)]);
print([outputDir 'Figure 4c XZCentersOfFaces.pdf'], '-dpdf');



%[Xgrid, Zgrid] = meshgrid(min(X_coordinate):fineGridResolution:max(X_coordinate), min(z_coordinate):fineGridResolution:max(z_coordinate));
[Xgrid, Zgrid] = meshgrid(linspace(min(X_coordinate), max(X_coordinate), 21), linspace(min(Z_coordinate), max(Z_coordinate), 20));
Ygrid = zeros(size(Xgrid)); % Set Y-coordinate to zero for all points


RBFq3TempG0_5 = griddata(X_coordinate, Z_coordinate, RBFsTempG0_5(:,3), Xgrid, Zgrid, 'nearest'); % Use cubic interpolation for smoother contours
RBFq3TempG1   = griddata(X_coordinate, Z_coordinate, RBFsTempG1(:,3), Xgrid, Zgrid, 'nearest'); % Use cubic interpolation for smoother contours
RBFq3TempG2   = griddata(X_coordinate, Z_coordinate, RBFsTempG2(:,3), Xgrid, Zgrid, 'nearest'); % Use cubic interpolation for smoother contours
RBFq3TempG2_5 = griddata(X_coordinate, Z_coordinate, RBFsTempG2_5(:,3), Xgrid, Zgrid, 'nearest'); % Use cubic interpolation for smoother contours

RBFq3TempM0_5 = griddata(X_coordinate, Z_coordinate, RBFsTempM0_5(:,3), Xgrid, Zgrid, 'nearest'); % Use cubic interpolation for smoother contours
RBFq3TempM1   = griddata(X_coordinate, Z_coordinate, RBFsTempM1(:,3), Xgrid, Zgrid, 'nearest'); % Use cubic interpolation for smoother contours
RBFq3TempM3   = griddata(X_coordinate, Z_coordinate, RBFsTempM3(:,3), Xgrid, Zgrid, 'nearest'); % Use cubic interpolation for smoother contours
RBFq3TempM7_5 = griddata(X_coordinate, Z_coordinate, RBFsTempM7_5(:,3), Xgrid, Zgrid, 'nearest'); % Use cubic interpolation for smoother contours


fontSize = 16; fontType = 'Times New Roman'; % Set the desired font type

figure(1)
TempG0_5 = surf(Xgrid, Zgrid, RBFq3TempG0_5, 'FaceAlpha',1);
grid on;  xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
zlabel('Gaussian RBF', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); % Gaussian 
eta_j_value = 0.5;
title(sprintf('\\eta_j = %.2f', eta_j_value), 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
xlim =[0,2]; ylim= [0,1.2];
colormap('jet');   % You can choose different colormaps
%shading interp;   % Interpolated shading for a smoother appearance
% Add a color bar (legend)
CC = colorbar; CC.FontSize = 14 ; CC.FontWeight = 'bold';
view(45, 30); % Adjust the view angle (optional)
dpi = 400; 
%print('GaussianTempG0_5.png', '-dpng', ['-r' num2str(dpi)]);
print([outputDir 'Figure 4a1 GaussianTempG0_5.png'], '-dpng', ['-r' num2str(dpi)]);

%print('GaussianTempG0_5.pdf', '-dpdf');

figure(2)
TempG1 = surf(Xgrid, Zgrid, RBFq3TempG1, 'FaceAlpha',1);
grid on;  xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
zlabel('Gaussian RBF', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); % Gaussian
eta_j_value = 1;
title(sprintf('\\eta_j = %.2f', eta_j_value), 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
xlim =[0,2]; ylim= [0,1.2];
colormap('jet');   % You can choose different colormaps
%shading interp;   % Interpolated shading for a smoother appearance
% Add a color bar (legend)
CC = colorbar; CC.FontSize = 10 ; CC.FontWeight = 'bold';
view(45, 30); % Adjust the view angle (optional)
dpi = 400; 
print([outputDir 'Figure 4a2 GaussianTempG1.png'], '-dpng', ['-r' num2str(dpi)]);
%print('GaussianTempG1.pdf', '-dpdf');

figure(3)
TempG2 = surf(Xgrid, Zgrid, RBFq3TempG2, 'FaceAlpha',1);
grid on;  xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
zlabel('Gaussian RBF', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); % Gaussian
eta_j_value = 2;
title(sprintf('\\eta_j = %.2f', eta_j_value), 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
xlim =[0,2]; ylim= [0,1.2];
colormap('jet');   % You can choose different colormaps
%shading interp;   % Interpolated shading for a smoother appearance
% Add a color bar (legend)
CC = colorbar; CC.FontSize = 10 ; CC.FontWeight = 'bold';
view(45, 30); % Adjust the view angle (optional)
dpi = 400; 
print([outputDir 'Figure 4a3 GaussianTempG2.png'], '-dpng', ['-r' num2str(dpi)]);
%print('GaussianTempG2.pdf', '-dpdf');

figure(4)
TempG2_5 = surf(Xgrid, Zgrid, RBFq3TempG2_5, 'FaceAlpha',1);
grid on;  xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
zlabel('Gaussian RBF', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); % Gaussian
eta_j_value = 2.5;
title(sprintf('\\eta_j = %.2f', eta_j_value), 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
xlim =[0,2]; ylim= [0,1.2];
colormap('jet');   % You can choose different colormaps
%shading interp;   % Interpolated shading for a smoother appearance
% Add a color bar (legend)
CC = colorbar; CC.FontSize = 10 ; CC.FontWeight = 'bold';
view(45, 30); % Adjust the view angle (optional)
dpi = 400; 
print([outputDir 'Figure 4a4 GaussianTempG2_5.png'], '-dpng', ['-r' num2str(dpi)]);
%print('GaussianTempG2_5.pdf', '-dpdf');

figure(5)
TempM0_5 = surf(Xgrid, Zgrid, RBFq3TempM0_5, 'FaceAlpha',1);
grid on;  xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
zlabel('Multiquadric RBF', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); % Multiquadric  
eta_j_value = 0.5;
title(sprintf('\\eta_j = %.2f', eta_j_value), 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
xlim =[0,2]; ylim= [0,1.2];
colormap('jet');   % You can choose different colormaps
%shading interp;   % Interpolated shading for a smoother appearance
% Add a color bar (legend)
CC = colorbar; CC.FontSize = 10 ; CC.FontWeight = 'bold';
view(45, 30); % Adjust the view angle (optional)
dpi = 400; 
print([outputDir 'Figure 4b1 MultiquadricTempM0_5.png'], '-dpng', ['-r' num2str(dpi)]);
%print('MultiquadricTempM0_5.pdf', '-dpdf');

figure(6)
TempM1 = surf(Xgrid, Zgrid, RBFq3TempM1, 'FaceAlpha',1);
grid on;  xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
zlabel('Multiquadric RBF', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); % Multiquadric  
eta_j_value = 1;
title(sprintf('\\eta_j = %.2f', eta_j_value), 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
xlim =[0,2]; ylim= [0,1.2];
colormap('jet');   % You can choose different colormaps
%shading interp;   % Interpolated shading for a smoother appearance
% Add a color bar (legend)
CC = colorbar; CC.FontSize = 10 ; CC.FontWeight = 'bold';
view(45, 30); % Adjust the view angle (optional)
dpi = 400; 
print([outputDir 'Figure 4b2 MultiquadricTempM1.png'], '-dpng', ['-r' num2str(dpi)]);
%print('MultiquadricTempM1.pdf', '-dpdf');

figure(7)
TempM3 = surf(Xgrid, Zgrid, RBFq3TempM3, 'FaceAlpha',1);
grid on;  xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
zlabel('Multiquadric RBF', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); % Multiquadric  
eta_j_value = 3;
title(sprintf('\\eta_j = %.2f', eta_j_value), 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
xlim =[0,2]; ylim= [0,1.2];
colormap('jet');   % You can choose different colormaps
%shading interp;   % Interpolated shading for a smoother appearance
% Add a color bar (legend)
CC = colorbar; CC.FontSize = 10 ; CC.FontWeight = 'bold';
view(45, 30); % Adjust the view angle (optional)
dpi = 400; 
print([outputDir 'Figure 4b3 MultiquadricTempM3.png'], '-dpng', ['-r' num2str(dpi)]);
%print('MultiquadricTempM3.pdf', '-dpdf');

figure(8)
TempM7_5 = surf(Xgrid, Zgrid, RBFq3TempM7_5, 'FaceAlpha',1);
grid on;  xlabel('X', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); ylabel('Z', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
zlabel('Multiquadric RBF', 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold'); % Multiquadric  
eta_j_value = 7.5;
title(sprintf('\\eta_j = %.2f', eta_j_value), 'FontSize', fontSize, 'FontName', fontType, 'FontWeight', 'bold');
xlim =[0,2]; ylim= [0,1.2];
colormap('jet');   % You can choose different colormaps
%shading interp;   % Interpolated shading for a smoother appearance
% Add a color bar (legend)
CC = colorbar; CC.FontSize = 10 ; CC.FontWeight = 'bold';
view(45, 30); % Adjust the view angle (optional)
dpi = 400; 
print([outputDir 'Figure 4b4 MultiquadricTempM7_5.png'], '-dpng', ['-r' num2str(dpi)]);
%print('MultiquadricTempM7_5.pdf', '-dpdf');



% Read images
img1 = imread([outputDir 'Figure 4a1 GaussianTempG0_5.png']);
img2 = imread([outputDir 'Figure 4a2 GaussianTempG1.png']);
img3 = imread([outputDir 'Figure 4a3 GaussianTempG2.png']);
img4 = imread([outputDir 'Figure 4a4 GaussianTempG2_5.png']);

img5 = imread([outputDir 'Figure 4b1 MultiquadricTempM0_5.png']);
img6 = imread([outputDir 'Figure 4b2 MultiquadricTempM1.png']);
img7 = imread([outputDir 'Figure 4b3 MultiquadricTempM3.png']);
img8 = imread([outputDir 'Figure 4b4 MultiquadricTempM7_5.png']);

% Combine first four images
combinedImg1GaussianRBF = [img1, img2; img3, img4]; % Adjust as needed

% Combine next four images
combinedImg2MultiquadricRBF = [img5, img6; img7, img8]; % Adjust as needed


% Save the first combined image
imwrite(combinedImg1GaussianRBF, [outputDir 'Figure 4a Combined1.png']);
figure;
imshow([outputDir 'Figure 4a Combined1.png']);
print([outputDir 'Figure 4a Combined1.pdf'], '-dpdf');

% Save the second combined image
imwrite(combinedImg2MultiquadricRBF, [outputDir 'Figure 4b Combined2.png']);

figure;
imshow([outputDir 'Figure 4b Combined2.png']);
print([outputDir 'Figure 4b Combined2.pdf'], '-dpdf');

