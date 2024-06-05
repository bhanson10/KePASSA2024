clear all; clc; close all;

%% Initialize Figure  
export = 1; % 0: don't export, 1: export
type   = 1; % 0: white figure, 1: black figure

if type==0
    createfigure_white; 
elseif type==1
    createfigure_black;
end

%% Particle Filter
p.type = 'sharp'; p.color = {[0, 1, 1], [0, 0.66, 0.66], [0, 0.33, 0.33]};  p.hgt = 1e-3; p.alpha = 1;  
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/PF/cmake-build-debug/Epochs/Europa/M1";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; 
for i=0:4:numFiles-1
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [x, P, n, t] = parse_nongaussian_txt(FILE_PATH);

    xest{count} = zeros(size(x(1,:)));
    for j=1:n
        xest{count} = xest{count}+x(j,:).*P(j);
    end

    for j=1:n
        x(j,:) = x(j,:) - xest{count};
    end

    X(:,1) = (count-1).*ones(size(x(:,1))); 
    X(:,2) = x(:,1);
    X(:,3) = x(:,2);
    plot_nongaussian_surface(X,P,[0.14,0.81,0.997],0.3,p);
    drawnow;

    count = count + 1; 
end
clear X; 

%% UKF
clear p; p.type = 'sharp'; p.color = {[1, 0, 0], [0.66, 0, 0], [0.33, 0, 0]};  p.hgt = 1e-3; p.alpha = 1;  
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/UKF/cmake-build-debug/Epochs/Europa/M0";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; 
for i=0:4:numFiles-1
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [~, Sig, t] = parse_gaussian_txt(FILE_PATH);

    X(1,1) = (count-1);
    X(2,1) = -0.03;
    X(3,1) = 0; 

    plot_gaussian_ellipsoid(X,Sig(1:2,1:2),[1,2,3],p);
    drawnow;

    count = count + 1;
end
clear X; 

%% GBEES
p.type = 'sharp'; p.color = {[0, 0.33, 0], [0, 0.66, 0], [0, 1, 0]}; p.hgt = 1e-3; p.alpha = 1;  
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/GBEES/cmake-build-debug/Epochs/Europa/M0";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; 
for i=0:4:numFiles-1
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [x, P, n, t] = parse_nongaussian_txt(FILE_PATH);

    X(:,1) = (count-1).*ones(size(x(:,1))); 
    X(:,2) = x(:,1);
    X(:,3) = x(:,2);
    plot_grid_surface(X,P,[0.011, 0.125, 0.58],0.3,p);
    drawnow;

    count = count + 1; clear X; 
end

if export
    if type==0
        export_fig('/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/png/Slides/title_slide_white_pcr3bp.png', figure(1), '-transparent', '-r800');
    elseif type==1
        export_fig('/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/png/Slides/title_slide_black_pcr3bp_2.png', figure(1), '-transparent', '-r800');
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
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
function [x, S, t] = parse_gaussian_txt(filename)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID)); 
    
    fgetl(fileID); % skip blank line

    line = split(fgetl(fileID)); % Read a line as a string
    x = [str2double(line{1});str2double(line{2});str2double(line{3});str2double(line{4})];
    
    fgetl(fileID); % skip blank line

    for i=1:4
        line = split(fgetl(fileID)); % Read a line as a string
        S(i,:) = [str2double(line{1});str2double(line{2});str2double(line{3});str2double(line{4})];
    end

    % Close the file
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createfigure_black

% Create figure
figure1 = figure('Color','none','Position',[300,1,900,865]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on'); grid(axes1,'on');

xlabel({'Time Past Epoch (hours)'},'FontName','Times','Rotation',-47, 'Position', [2.248488234361811,-0.056402222672616,-0.00270023732813]);

view(axes1,[70 20]);

xlim([0,4.1]);

set(axes1,'Color','none','DataAspectRatio',[140 5 1],'FontName','Times',...
    'FontSize',26,'GridColor',[0 0 0],'LineWidth',2,'XColor',[1 1 1],'XTick',...
    [0 1 2 3 4.1],'XTickLabel',{'0','3.5','7','10.5','14'},'YColor',[1 1 1],...
    'YTick',[],'ZColor',[1 1 1],'ZTick',[]);
drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createfigure_white

% Create figure
figure1 = figure('Position',[300,1,900,865]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on'); grid(axes1,'on');

xlabel({'Time Past Epoch (hours)'},'FontName','Times','Rotation',-47, 'Position', [2.248488234361811,-0.056402222672616,-0.00270023732813]);

view(axes1,[70 20]);

xlim([0,4.1]);

set(axes1,'DataAspectRatio',[140 5 1],'FontName','Times',...
    'FontSize',26,'LineWidth',2,'XTick',...
    [0 1 2 3 4.1],'XTickLabel',{'0','3.5','7','10.5','14'},...
    'YTick',[],'ZColor',[1 1 1],'ZTick',[]);
drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%