clear all; clc; close all; 

%% Initial Condition
const.LU = 668519; const.TU = 48562; const.cV = [const.LU, const.LU, const.LU/const.TU, const.LU/const.TU]; 
lbl.XString = 'Time Past Epoch (hours)'; lbl.YString = 'Resolution';
initialize_figures('n', 1, 'lgd', 1, 'lbl', lbl, 'spacing', {[183,310,1083,542]});

%% Particle Filter
p.plt = 0;
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/PF/cmake-build-debug/Epochs/Europa/M1";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; 
for i=[0:numFiles-1]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [x_pf, P_pf, n_pf, t_pf(count)] = parse_nongaussian_txt(FILE_PATH);

    % Converting units
    x_pf(:,1:4) = x_pf(:,1:4).*const.cV; 
    t_pf(count) = t_pf(count)*const.TU;

    [~,shp_p] = plot_nongaussian_surface(x_pf(:,1:2),P_pf,[0.99925],0.3,p);

    resolution_pf(count) = n_pf/area(shp_p{1}); 
    
    count = count + 1;
end

figure(1); 
plot(t_pf/3600, resolution_pf, '-o', 'Color', [0, 0.66, 0.66], 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'PF')

%% GBEES
p.plt = 0;
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/GBEES/cmake-build-debug/Epochs/Europa/M0";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1;
for i=[0:numFiles-1]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(FILE_PATH);

    % Converting units
    x_gbees(:,1:4) = x_gbees(:,1:4).*const.cV; 
    t_gbees(count) = t_gbees(count)*const.TU;

    [~,shp_p] = plot_grid_surface(x_gbees(:,1:2),P_gbees,[0.011, 0.125, 0.58],0.3,p);

    resolution_gbees(count) = n_gbees/area(shp_p{1}); 
    
    count = count + 1;
end

figure(1); 
plot(t_gbees/3600, resolution_gbees, '-square', 'Color', 'g', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'GBEES')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        x(count, :) = [str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5})];
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%