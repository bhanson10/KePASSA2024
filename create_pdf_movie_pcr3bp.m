clear all; clc; close all;

%% Initialize Figure  
export = 1; % 0: don't export, 1: export
type   = 0; % 0: white figure, 1: black figure

if type==0
    createfigure_white; 
elseif type==1
    createfigure_black;
end

%% Particle Filter
p.type = 'sharp'; p.color = {[0, 1, 1], [0, 0.66, 0.66], [0, 0.33, 0.33]};  p.hgt = 1e-3; p.alpha = 1;  
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/PF/cmake-build-debug/Epochs/Europa/M0";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; 
for i=0:numFiles-1
     if i~=0
        for j=1:numel(plots)
            delete(plots{j});
        end
     end

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
    
    p.color = {[0, 1, 1], [0, 0.66, 0.66], [0, 0.33, 0.33]}; p.alpha = 1; 
    [plots, ~] = plot_nongaussian_surface(X(:,1:3),P,[0.14,0.81,0.997],0.5,p);
    p.color = [0, 1, 1]; p.alpha = 0.05;
    plot_nongaussian_surface(X(:,1:3),P,[0.997],0.5,p);

    drawnow;

    frames(count) = getframe(gcf); 

    count = count + 1; 
end

if export
    create_video(frames,'/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/animations/europa_pf_pdf_movie.mp4',20);
end

%% UKF
%{
p.color = 'green'; 
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/PCR3BP/UKF/cmake-build-debug/Epochs/Europa/M0";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; 
for i=0:numFiles-1
    if i~=0
        for j=1:numel(plots)
            delete(plots{j});
        end
    end

    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [~, Sig, t] = parse_gaussian_txt(FILE_PATH);

    figure(1);
    X(1,1) = (count-1)*5.95;
    X(2,1) = 0;
    X(3,1) = 0;
    S = Sig(1:3,1:3);

    p.alpha = [0.8, 0.4, 0.25]; 
    [plots, shps] = plot_gaussian_ellipsoid(X,S([3,2,1],[3,2,1]),[1,2,3],p);
    p.alpha = 0.05; 
    plot_gaussian_ellipsoid(X,S([3,2,1],[3,2,1]),[3],p);
    drawnow;

    frames(count) = getframe(gcf); 

    count = count + 1; 
end

if export
    create_video(frames,'/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/animations/europa_ukf_pdf_movie_4.mp4',20);
end
%}
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

set(axes1,'Color','none','DataAspectRatio',[1 1 1],'FontName','Times',...
    'FontSize',26,'GridColor',[0 0 0],'LineWidth',2,'XColor',[1 1 1],'XTick',...
    [0 1 2 3 4.1],'XTickLabel',{'0','3.5','7','10.5','14'},'YColor',[1 1 1],...
    'YTick',[],'ZColor',[1 1 1],'ZTick',[]);
drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createfigure_white

% Create figure
figure1 = figure('Color','white', 'Position',[300,1,900,865]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on'); grid(axes1,'on');

xlabel({'Time Past Epoch (hours)'},'FontName','Times','Rotation',-47, 'Position', [2.248488234361811,-0.056402222672616,-0.00270023732813]);

view(axes1,[70 20]);

%xlim([0,4.1]);

set(axes1,'DataAspectRatio',[1 1 1],'FontName','Times',...
    'FontSize',26,'LineWidth',2,'XTick',...
    [0 1 2 3 4.1],'XTickLabel',{'0','3.5','7','10.5','14'},...
    'YTick',[],'ZColor',[1 1 1],'ZTick',[]);
drawnow; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%