clear all; clc; close all;

%% Initialize Figure  
export = 0; % 0: don't export, 1: export
type   = 0; % 0: white figure, 1: black figure

if type==0
    createfigure_white; 
elseif type==1
    createfigure_black;
end

%% Particle Filter
p.type = 'sharp'; p.color = 'cyan'; p.alpha = [0.8, 0.4, 0.25]; 
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/PF/cmake-build-debug/Epochs/Europa/Classical/Revs6";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; 
for i=0:4:16
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 
    fileID = fopen(FILE_PATH, 'r'); 
    t = str2double(fgetl(fileID)); 

    [x, P, n, t] = parse_nongaussian_txt(FILE_PATH, 'coe');

    xest{count} = zeros(size(x(1,:)));
    for j=1:n
        xest{count} = xest{count}+x(j,:).*P(j);
    end

    for j=1:n
        x(j,:) = x(j,:) - xest{count};
    end

    %Scale the data so it is square
    xrange = max(x(:,1)) - min(x(:,1));
    yrange = max(x(:,2)) - min(x(:,2));
    zrange = max(x(:,3)) - min(x(:,3));

    max_range = max([xrange, yrange, zrange]);

    xscale{count} = max_range/xrange;
    yscale{count} = max_range/yrange;
    zscale{count} = 1.15*max_range/zrange;

    figure(1); 
    X(:,1) = x(:,3).*zscale{count} + (count-1)*800;
    X(:,2) = x(:,2).*yscale{count};
    X(:,3) = x(:,1).*xscale{count};

    plot_nongaussian_surface(X(:,1:3),P,[0.14,0.81,0.997],1,p);
    drawnow;

    count = count + 1; fclose(fileID);

    if i==0
        if type==0
            figure2 = figure('Color','white'); hold on; grid off; view(45,45); axis square;  
            set(gca,'Color','white');
        elseif type==1
            figure2 = figure('Color','none'); hold on; grid off; view(45,45); axis square;  
            set(gca,'Color','none');
        end
        light('Parent',gca,'Position',[1 0 0]);
        plot_nongaussian_surface(X(:,1:3),P,[0.14,0.81,0.997],0.05,p);
    end
end
clear X; 

%% UKF 
clear p; p.color = 'red'; p.alpha = [0.8, 0.4, 0.25]; 
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/UKF/cmake-build-debug/Epochs/Europa/Equinoctial/Revs6";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; 
for i=0:4:16
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 
    fileID = fopen(FILE_PATH, 'r'); 
    t = str2double(fgetl(fileID)); 

    [~, Sig, t] = parse_gaussian_txt(FILE_PATH, 'eoe');

    figure(1);
    X(1,1) = (count-1)*800;
    X(2,1) = -2000;
    X(3,1) = 0;

    scale = [xscale{count}, yscale{count}, zscale{count}];
    S = Sig(1:3,1:3).*(scale'*scale);  

    plot_gaussian_ellipsoid(X,S([3,2,1],[3,2,1]),[1,2,3],p);
    drawnow;

    count = count + 1; fclose(fileID);

    if i==0
        if type==0
            figure3 = figure('Color','white'); hold on; grid off; view(45,45); axis square;  
            set(gca,'Color','white');
        elseif type==1
            figure3 = figure('Color','none'); hold on; grid off; view(45,45); axis square;  
            set(gca,'Color','none');
        end
        light('Parent',gca,'Position',[1 0 0]);
        plot_gaussian_ellipsoid(X,S([3,2,1],[3,2,1]),[1,2,3],p);
    end
end

if export
    if type==0
        export_fig('/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/png/Slides/title_slide_white_europa_r2bp.png', figure(1), '-transparent', '-r800');
    elseif type==1
        export_fig('/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/png/Slides/title_slide_black_europa_r2bp.png', figure(1), '-transparent', '-r800');
    end
    export_fig('/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/png/Slides/title_slide_pf_europa_r2bp.png', figure(2), '-transparent', '-r800');
    export_fig('/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/png/Slides/title_slide_ekf_europa_r2bp.png', figure(3), '-transparent', '-r800');
end
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
figure1 = figure('Color','none', 'Position',[169,88,1206,694]); 

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create light
light('Parent',axes1,'Position',[1 0 0]);

% Create xlabel
xlabel('Revolutions','FontSize',30,'FontName','Times','Rotation',-6.5,...
    'Interpreter','latex', 'Position',[1300,2000,-3500]);

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0,Inf]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-Inf Inf]);
% Uncomment the following line to preserve the Z-limits of the axes
zlim(axes1,[-600 300]);
view(axes1,[13.134728311494872,28.439219287760807]);
grid(axes1,'on');

% Set the remaining axes properties
set(axes1,'Color','none','DataAspectRatio',[1 1 1],'FontName','Times',...
    'FontSize',26,'GridColor',[0 0 0],'LineWidth',2,'XColor',[1 1 1],'XTick',...
    [0 800 1600 2400 3200],'XTickLabel',{'0','1','2','3','4'},...
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
    'Interpreter','latex', 'Position',[1300,2000,-3500]);

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0,Inf]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[-Inf Inf]);
% Uncomment the following line to preserve the Z-limits of the axes
zlim(axes1,[-600 300]);
view(axes1,[13.134728311494872,28.439219287760807]);
grid(axes1,'on');

% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'FontName','Times',...
    'FontSize',26,'LineWidth',2,'XTick',...
    [0 800 1600 2400 3200],'XTickLabel',{'0','1','2','3','4'},...
    'YTick',[],'ZTick',[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%