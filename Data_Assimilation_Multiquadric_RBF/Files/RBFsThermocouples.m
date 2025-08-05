%################## Please read before running

% Reading the coordinates of thermocouples through thermocouplesDict. 
% For defining the location of RBFs, we project 5/9 or more selected thermocouples among 25/50/75/100 of them on the hotSide surface. 
% The center of each function is at the projection of these selected thermocouples on the boundary hotSide
% RBFs must be less than measurements or at least are equal to measurements.

% Plotting the location of all thermocouples and visually selecting some of them using for projection on the hotside as RBFs
% Adding a plot of the cube with transparent surfaces,

%################## Please read before running

clc; clearvars; close all;

% Set the current directory to "my_folder" on the desktop


% Define the path to the file
filePath = './constant/thermocouplesDict';

% Read the file
fid = fopen(filePath, 'r');
if fid == -1
    error('Unable to open the file.');
end

% Initialize arrays to store coordinates
x = [];
y = [];
z = [];

% Read the lines of the file
while ~feof(fid)
    line = fgetl(fid);
    
    % Check if the line contains coordinates
    if contains(line, '(') && contains(line, ')')
        % Extract coordinates from the line
        C = textscan(line, '(%f %f %f)');
        x = [x; C{1}];
        y = [y; C{2}];
        z = [z; C{3}];
    end
end
% Close the file
fclose(fid);

% Create the first figure and plot all thermocouples
fig1 = figure;
scatter3(x, y, z, 'filled'); 
%xlabel('X'); ylabel('Y'); zlabel('Z');title('All Thermocouples');grid on;

% Adding a plot of the cube with transparent surfaces, et the FaceAlpha property to make the surfaces transparent.
% Define the dimensions of the cube
cubeX = [0, 2, 2, 0, 0, 2, 2, 0];
cubeY = [0, 0, 0.1, 0.1, 0, 0, 0.1, 0.1];
cubeZ = [0, 0, 0, 0, 1.2, 1.2, 1.2, 1.2];
 
% Create a patch for the cube and set transparency
FaceAlphaTransparency = 0.04;
patch([cubeX(1), cubeX(2), cubeX(6), cubeX(5)], [cubeY(1), cubeY(2), cubeY(6), cubeY(5)], [cubeZ(1), cubeZ(2), cubeZ(6), cubeZ(5)], 'b', 'FaceAlpha', FaceAlphaTransparency);
hold on;
patch([cubeX(2), cubeX(3), cubeX(7), cubeX(6)], [cubeY(2), cubeY(3), cubeY(7), cubeY(6)], [cubeZ(2), cubeZ(3), cubeZ(7), cubeZ(6)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(3), cubeX(4), cubeX(8), cubeX(7)], [cubeY(3), cubeY(4), cubeY(8), cubeY(7)], [cubeZ(3), cubeZ(4), cubeZ(8), cubeZ(7)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(4), cubeX(1), cubeX(5), cubeX(8)], [cubeY(4), cubeY(1), cubeY(5), cubeY(8)], [cubeZ(4), cubeZ(1), cubeZ(5), cubeZ(8)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(1), cubeX(2), cubeX(3), cubeX(4)], [cubeY(1), cubeY(2), cubeY(3), cubeY(4)], [cubeZ(1), cubeZ(2), cubeZ(3), cubeZ(4)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(5), cubeX(6), cubeX(7), cubeX(8)], [cubeY(5), cubeY(6), cubeY(7), cubeY(8)], [cubeZ(5), cubeZ(6), cubeZ(7), cubeZ(8)], 'b', 'FaceAlpha', FaceAlphaTransparency);
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Visualization of a Transparent Mold with Thermocouple Locations (Blue Dots)');
grid on;

% Set the DataAspectRatio to manual to plot in real-world coordinates
axis equal;
close(fig1);
%%
% Find the indices of the corner and middle thermocouples
%cornerMiddleIndices = [1, 5, 13, 21, 25];  % 5 out of 25
cornerMiddleIndices = [1, 10, 55, 91, 100];   % 5 out of 100

% Create the second figure and plot the five corner and middle thermocouples in red and the others in blue
fig2 = figure;
% Maximize the figure window
set(fig2, 'WindowState', 'maximized');  % Works in newer MATLAB versions (R2020b+)

scatter3(x(cornerMiddleIndices), zeros(length(cornerMiddleIndices)), z(cornerMiddleIndices), 'filled', 'r');
hold on;
otherIndices = setdiff(1:numel(x), cornerMiddleIndices);
scatter3(x(otherIndices), y(otherIndices), z(otherIndices), 'filled', 'b');
hold off;xlabel('X');ylabel('Y');zlabel('Z'); title('Visualization of a Transparent Mold with Thermocouple Locations (Blue Dots) and RBF Centers (Red Dots)');grid on;

patch([cubeX(1), cubeX(2), cubeX(6), cubeX(5)], [cubeY(1), cubeY(2), cubeY(6), cubeY(5)], [cubeZ(1), cubeZ(2), cubeZ(6), cubeZ(5)], 'b', 'FaceAlpha', FaceAlphaTransparency);
hold on;
patch([cubeX(2), cubeX(3), cubeX(7), cubeX(6)], [cubeY(2), cubeY(3), cubeY(7), cubeY(6)], [cubeZ(2), cubeZ(3), cubeZ(7), cubeZ(6)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(3), cubeX(4), cubeX(8), cubeX(7)], [cubeY(3), cubeY(4), cubeY(8), cubeY(7)], [cubeZ(3), cubeZ(4), cubeZ(8), cubeZ(7)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(4), cubeX(1), cubeX(5), cubeX(8)], [cubeY(4), cubeY(1), cubeY(5), cubeY(8)], [cubeZ(4), cubeZ(1), cubeZ(5), cubeZ(8)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(1), cubeX(2), cubeX(3), cubeX(4)], [cubeY(1), cubeY(2), cubeY(3), cubeY(4)], [cubeZ(1), cubeZ(2), cubeZ(3), cubeZ(4)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(5), cubeX(6), cubeX(7), cubeX(8)], [cubeY(5), cubeY(6), cubeY(7), cubeY(8)], [cubeZ(5), cubeZ(6), cubeZ(7), cubeZ(8)], 'b', 'FaceAlpha', FaceAlphaTransparency);

% Set the DataAspectRatio to manual to plot in real-world coordinates
axis equal;


% Tighten layout manually (alternative to 'Bounds', 'tight')
set(gca, 'LooseInset', get(gca, 'TightInset'));
set(gcf, 'Color', 'none');  % Optional: remove figure background

% Define path and save with high resolution
outputPath = fullfile('..', 'Results', 'Figure 2 Thermocouples_RBF_Centers.png');
print(gcf, outputPath, '-dpng', '-r300');  % -dpng = PNG format, -r300 = 300 DPI


% Find the indices of corner, middle and diagonal thermocouples

cornerMiddleDiagonalIndices = [1, 10, 12, 19, 23, 28, 34, 37, 45, 46, 55, 56, 64, 66, 73, 78, 82, 89, 91, 100];   % 20 out of 100
%cornerMiddleDiagonalIndices = [1, 5, 7, 9, 13, 17, 19, 21, 25]; % 9 out of 25
% Create the third figure with nine corner, middle and diagonal thermocouples in red and the others in blue
fig3= figure;
scatter3(x(cornerMiddleDiagonalIndices), zeros(length(cornerMiddleDiagonalIndices)), z(cornerMiddleDiagonalIndices), 'filled', 'r');
hold on;
otherIndices1 = setdiff(1:numel(x), cornerMiddleDiagonalIndices);
scatter3(x(otherIndices1), y(otherIndices1), z(otherIndices1), 'filled', 'b');
hold off; xlabel('X'); ylabel('Y'); zlabel('Z'); title('Visualization of a Transparent Mold with Thermocouple Locations (Blue Dots) and RBF Centers (Red Dots)'); grid on;patch([cubeX(1), cubeX(2), cubeX(6), cubeX(5)], [cubeY(1), cubeY(2), cubeY(6), cubeY(5)], [cubeZ(1), cubeZ(2), cubeZ(6), cubeZ(5)], 'b', 'FaceAlpha', FaceAlphaTransparency);
hold on;
patch([cubeX(2), cubeX(3), cubeX(7), cubeX(6)], [cubeY(2), cubeY(3), cubeY(7), cubeY(6)], [cubeZ(2), cubeZ(3), cubeZ(7), cubeZ(6)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(3), cubeX(4), cubeX(8), cubeX(7)], [cubeY(3), cubeY(4), cubeY(8), cubeY(7)], [cubeZ(3), cubeZ(4), cubeZ(8), cubeZ(7)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(4), cubeX(1), cubeX(5), cubeX(8)], [cubeY(4), cubeY(1), cubeY(5), cubeY(8)], [cubeZ(4), cubeZ(1), cubeZ(5), cubeZ(8)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(1), cubeX(2), cubeX(3), cubeX(4)], [cubeY(1), cubeY(2), cubeY(3), cubeY(4)], [cubeZ(1), cubeZ(2), cubeZ(3), cubeZ(4)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(5), cubeX(6), cubeX(7), cubeX(8)], [cubeY(5), cubeY(6), cubeY(7), cubeY(8)], [cubeZ(5), cubeZ(6), cubeZ(7), cubeZ(8)], 'b', 'FaceAlpha', FaceAlphaTransparency);

% Set the DataAspectRatio to manual to plot in real-world coordinates
axis equal;
close(fig3);

% Find the indices of the corner and middle thermocouples
% cornerMiddleIndices = [3, 12, 14, 23];  % 4 out of 25
cornerMiddleIndices = [13, 18, 43, 48, 83, 88];  % 6 out of 100
% Create the second figure and plot the five corner and middle thermocouples in red and the others in blue
fig4=figure;
scatter3(x(cornerMiddleIndices), zeros(length(cornerMiddleIndices)), z(cornerMiddleIndices), 'filled', 'r');
hold on;
otherIndices = setdiff(1:numel(x), cornerMiddleIndices);
scatter3(x(otherIndices), y(otherIndices), z(otherIndices), 'filled', 'b');
hold off;
xlabel('X', 'FontName','Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Y', 'FontName','Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('Z', 'FontName','Times New Roman', 'FontSize', 14, 'FontWeight', 'bold'); 
%title('Visualization of a Transparent Mold with Thermocouple Locations (Blue Dots) and RBF Centers (Red Dots)');grid on;

patch([cubeX(1), cubeX(2), cubeX(6), cubeX(5)], [cubeY(1), cubeY(2), cubeY(6), cubeY(5)], [cubeZ(1), cubeZ(2), cubeZ(6), cubeZ(5)], 'b', 'FaceAlpha', FaceAlphaTransparency);
hold on;
patch([cubeX(2), cubeX(3), cubeX(7), cubeX(6)], [cubeY(2), cubeY(3), cubeY(7), cubeY(6)], [cubeZ(2), cubeZ(3), cubeZ(7), cubeZ(6)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(3), cubeX(4), cubeX(8), cubeX(7)], [cubeY(3), cubeY(4), cubeY(8), cubeY(7)], [cubeZ(3), cubeZ(4), cubeZ(8), cubeZ(7)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(4), cubeX(1), cubeX(5), cubeX(8)], [cubeY(4), cubeY(1), cubeY(5), cubeY(8)], [cubeZ(4), cubeZ(1), cubeZ(5), cubeZ(8)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(1), cubeX(2), cubeX(3), cubeX(4)], [cubeY(1), cubeY(2), cubeY(3), cubeY(4)], [cubeZ(1), cubeZ(2), cubeZ(3), cubeZ(4)], 'b', 'FaceAlpha', FaceAlphaTransparency);
patch([cubeX(5), cubeX(6), cubeX(7), cubeX(8)], [cubeY(5), cubeY(6), cubeY(7), cubeY(8)], [cubeZ(5), cubeZ(6), cubeZ(7), cubeZ(8)], 'b', 'FaceAlpha', FaceAlphaTransparency);

% Set the DataAspectRatio to manual to plot in real-world coordinates
axis equal;
close(fig4);