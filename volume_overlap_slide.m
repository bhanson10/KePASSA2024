clear all; clc; close all; 

% Example: 
%% Initial Condition
rng(2024);
traj.sys = 'R2BP'; name = "Earth"; nm = 0;
filename = append('./ICs/', traj.sys,'/IC_', name,'.csv');
po = readmatrix(filename); const.mu = po(8);
initialize_figures('n', 1:4, 'bg', 'clear', 'lght', 1, 'vw', {[45,45]}, 'spacing', {[251,84,903,782]});

%% Particle Filter
p.type = 'sharp'; p.color = 'blue'; p.means = 0; p.display = 0; p.name = "PF"; p.alpha = [0.1]; 
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/PF/cmake-build-debug/Epochs/Moon/Classical/Revs4";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);
FILE_PATH = DATA_PATH + "/pdf_" + num2str(numFiles-1) + ".txt"; 
[x_pf, P_pf, n_pf, t_pf] = parse_nongaussian_txt(FILE_PATH, 'coe');
figure(1);
plot_nongaussian_surface(x_pf(:,1:3),P_pf,[0.997],0.25,p);
figure(2);
[~,shp_pf] =  plot_nongaussian_surface(x_pf(:,1:3),P_pf,[0.997],0.25,p);
shp1 = shp_pf{1};

% Metric Slide
figure(3);
plot_nongaussian_surface(x_pf(:,1:3),P_pf,[0.2,0.74,0.971],0.25,p);
colors = jet(100);
scaled_w = (P_pf - min(P_pf))./(max(P_pf)-min(P_pf));
colors = colors(round(scaled_w.*99)+1,:);
figure(4)
n = 1000; 
scatter3(x_pf(1:n,1), x_pf(1:n,2), x_pf(1:n,3), 20, colors(1:n,:),'filled',"HandleVisibility","off");


%% UKF
p.type = 'sharp'; p.color = 'red'; p.means = 0; p.display = 0; p.name = "PF"; p.alpha = [0.2]; 
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/UKF/cmake-build-debug/Epochs/Moon/Classical/Revs4";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);
FILE_PATH = DATA_PATH + "/pdf_" + num2str(numFiles-1) + ".txt"; 
[x_ukf, S_ukf, t_ukf] = parse_gaussian_txt(FILE_PATH, 'coe');
sd = 3; 
figure(1)
[~,shp_ukf] = plot_gaussian_ellipsoid(x_ukf(1:3),S_ukf(1:3,1:3),sd,p);
figure(2)
plot_gaussian_ellipsoid(x_ukf(1:3),S_ukf(1:3,1:3),sd,p);
shp2.Points = shp_ukf{1};
shp2.mean   = x_ukf(1:3); 
shp2.cov    = S_ukf(1:3,1:3); 
shp2.sd     = sd;

isAlpha1 = isa(shp1, 'alphaShape'); 
isAlpha2 = isa(shp2, 'alphaShape'); 
if ~exist('npts', 'var'), npts = 40; end

x_min = min(shp1.Points(:,1)); x_max = max(shp1.Points(:,1));
y_min = min(shp1.Points(:,2)); y_max = max(shp1.Points(:,2));
z_min = min(shp1.Points(:,3)); z_max = max(shp1.Points(:,3));

x_min = min(x_min,min(shp2.Points(:,1))); x_max = max(x_max,max(shp2.Points(:,1)));
y_min = min(y_min,min(shp2.Points(:,2))); y_max = max(y_max,max(shp2.Points(:,2)));
z_min = min(z_min,min(shp2.Points(:,3))); z_max = max(z_max,max(shp2.Points(:,3)));

count = 0; count_union = 0; count_intersection = 0; count_in_1 = 0; count_in_2 = 0; 
for i=linspace(x_min,x_max,npts)
    for j=linspace(y_min,y_max,npts)
        for k=linspace(z_min,z_max,npts)
            pt = [i, j, k]; count = count + 1; 
            pts(count, :) = pt;

            % Checking in Shp1
            inShp1 = 0; 
            if isAlpha1
                if inShape(shp1,pt) 
                    inShp1 = 1; count_in_1 = count_in_1 + 1; 
                    pts_in_1(count_in_1,:) = pt; 
                end
            else
                dM = sqrt((pt' - shp1.mean)'*inv(shp1.cov)*(pt' - shp1.mean));
                if dM <= shp1.sd
                    inShp1 = 1; count_in_2 = count_in_1 + 1; 
                    pts_in_1(count_in_1,:) = pt; 
                end
            end

            % Checking in Shp2
            inShp2 = 0; 
            if isAlpha2
                if inShape(shp2,pt) 
                    inShp2 = 1; count_in_2 = count_in_2 + 1; 
                    pts_in_2(count_in_2,:) = pt; 
                end
            else
                dM = sqrt((pt' - shp2.mean)'*inv(shp2.cov)*(pt' - shp2.mean));
                if dM <= shp2.sd
                    inShp2 = 1; count_in_2 = count_in_2 + 1; 
                    pts_in_2(count_in_2,:) = pt;  
                end
            end

            if inShp1||inShp2
                count_union = count_union + 1; 
            end

            if inShp1&&inShp2
                count_intersection = count_intersection + 1; 
                pts_intersection(count_intersection, :) = pt;
            end
        end
    end
end


%scatter3(pts(:,1),pts(:,2),pts(:,3),5,'filled','k');
scatter3(pts_in_1(:,1),pts_in_1(:,2),pts_in_1(:,3),30,'filled','b');
scatter3(pts_in_2(:,1),pts_in_2(:,2),pts_in_2(:,3),30,'filled','r');
scatter3(pts_intersection(:,1),pts_intersection(:,2),pts_intersection(:,3),60,'filled','magenta');
j = count_intersection / count_union,

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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