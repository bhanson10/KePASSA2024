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
p.type = 'sharp'; p.color = 'cyan';  
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/PF2/cmake-build-debug/Epochs/Europa/Classical/Revs_movie4";
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

    [x, P, n, t] = parse_nongaussian_txt(FILE_PATH, 'coe');

    xest{count} = zeros(size(x(1,:)));
    for j=1:n
        xest{count} = xest{count}+x(j,:).*P(j);
    end

    for j=1:n
        x(j,:) = x(j,:) - xest{count};
    end

    X(:,1) = x(:,3) + (count-1)*5.95; % the spacing should be enough such that the figure is the same size as the title one
    X(:,2) = x(:,2);
    X(:,3) = x(:,1);

    p.alpha = [0.8, 0.4, 0.25]; 
    [plots, ~] = plot_nongaussian_surface(X(:,1:3),P,[0.14,0.81,0.997],0.5,p);
    p.alpha = 0.05;
    plot_nongaussian_surface(X(:,1:3),P,[0.997],0.5,p);

    drawnow;

    frames(count) = getframe(gcf); 

    count = count + 1; 
end

if export
    create_video(frames,'/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/animations/europa_pf_pdf_movie.mp4',20);
end

%% EKF
%{
p.color = 'green'; 
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/UKF2/cmake-build-debug/Epochs/Europa/Equinoctial/Revs_movie4";
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

    [~, Sig, t] = parse_gaussian_txt(FILE_PATH, 'eoe');

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
function [x, P, n, t] = parse_nongaussian_txt(filename,oetype)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        oe = [str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5});str2double(line{6});mod(str2double(line{7}),2*pi)];
        x(count, :) = oe2rv(oe,398600.4354,'rad','M',oetype); % Convert to Cartesian
        count = count + 1; 
    end
    
    % Close the file
    fclose(fileID);
    n = length(P); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, S, t] = parse_gaussian_txt(filename, oetype)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID)); 
    
    fgetl(fileID); % skip blank line

    line = split(fgetl(fileID)); % Read a line as a string
    oe = [str2double(line{1});str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5});str2double(line{6})];

    if ~strcmp(oetype, 'cart')
        x  = oe2rv(oe,398600.4354,'rad','M',oetype); % Convert to Cartesian
    else
        x = oe; 
    end

    fgetl(fileID); % skip blank line

    for i=1:6
        line = split(fgetl(fileID)); % Read a line as a string
        S(i,:) = [str2double(line{1});str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5});str2double(line{6})];
    end
    
    if ~strcmp(oetype, 'cart')
        S = oe2rv_cov(oe,S,398600.4354,'M',oetype);
    end

    % Close the file
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createfigure_black
figure1 = figure('Color','black','Position',[169,88,1206,694]); 

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create light
light('Parent',axes1,'Position',[1 0 0]);

% Create xlabel
xlabel('Revolutions','FontSize',30,'FontName','Times','Rotation',-6.5,...
    'Interpreter','latex');

xlim(axes1,[-3,inf]);
% ylim(axes1,[-Inf Inf]);
% zlim(axes1,[-600 300]);
view(axes1,[13.134728311494872,28.439219287760807]);
grid(axes1,'on');

% Set the remaining axes properties
set(axes1,'Color','black','DataAspectRatio',[1 1 1],'FontName','Times',...
    'FontSize',26,'GridColor',[0 0 0],'LineWidth',2,'XColor',[1 1 1],'XTick',...
    [0 600 1200 1800 2400 3000 3600],'XTickLabel',{'0','1','2','3','4','5','6'},...
    'YColor',[1 1 1],'YTick',[],'ZColor',[1 1 1],'ZTick',[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function createfigure_white
figure1 = figure('Position',[169,88,1206,694]); 

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create light
light('Parent',axes1,'Position',[1 0 0]);

% Create xlabel
xlabel('Revolutions','FontSize',30,'FontName','Times','Rotation',-6.5,...
    'Interpreter','latex');

xlim(axes1,[-3,inf]);
% ylim(axes1,[-Inf Inf]);
% zlim(axes1,[-600 300]);
view(axes1,[13.134728311494872,28.439219287760807]);
grid(axes1,'on');

count = 1; 
for i=0:0.25:4
    ticks(count) = i*1190; 
    ticklabels{count} = num2str(i); 
    count = count + 1; 
end

% Set the remaining axes properties
set(axes1,'Color','white','DataAspectRatio',[1 1 1],'FontName','Times',...
    'FontSize',26,'GridColor',[0 0 0],'LineWidth',2,'XColor',[0 0 0],'XTick',...
    ticks,'XTickLabel',ticklabels,...
    'YColor',[0 0 0],'YTick',[],'ZColor',[0 0 0],'ZTick',[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%