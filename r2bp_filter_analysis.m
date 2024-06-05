clear all; close all; clc;

% r2bp_filter_analysis.m
% Benjamin Hanson, 2024

%% Initial Condition
rng(2024);
traj.sys = 'R2BP'; name = 'Europa_2';
filename = append('./ICs/', traj.sys,'/IC_', name,'.csv');
po = readmatrix(filename); const.mu = po(8);
traj.rv0=[po(1); po(2); po(3); po(4); po(5); po(6)]; 
traj.oe0=rv2oe(traj.rv0, const.mu, 'rad', 'M');
const.T = po(7); const.r = po(9); TOL = 1e-5; 
colors = distinguishable_colors(6,{'w','y'});

%% Initializing Figures
lbl.XString = '$x$ (km)';
lbl.YString = '$y$ (km)';
lbl.ZString = '$z$ (km)';
initialize_figures('n', 1, 'spacing', {[50 100 700 700]},'lbl', lbl, 'vw', {[65,10]}, 'axs', {'equal'}); 

% Loading Europa
imageFile = '/Users/bhanson/OneDrive - UC San Diego/UCSD/Conferences/KePASSA2024/Figures/png/Slides/Europa_Map.png';
europaImg = imread(imageFile);

[X, Y, Z] = sphere(100);  
X = X.*const.r; Y = Y.*const.r; Z = Z.*const.r;
surf(X, Y, Z, 'FaceColor', 'texturemap', 'CData', europaImg, 'EdgeColor', 'none', 'HandleVisibility', 'off');
drawnow;

lbl.XString = '$v_x$ (km/s)';
lbl.YString = '$v_y$ (km/s)';
lbl.ZString = '$v_z$ (km/s)';
initialize_figures('n', 2, 'spacing', {[800 100 700 700]}, 'lbl', lbl, 'vw', {[55,10]}, 'axs', {'equal'}); 

%% Truth
rev = 4; 
num_dist = 4; 
measure_t = const.T/num_dist; 
tspan = [0:measure_t:rev*const.T]; 
t = 0; x0 = traj.oe0; x = x0'; t_idx = 1;
options = odeset('MaxStep', 5, 'InitialStep', 0.1, 'RelTol', 1e-6);
count = 1; 
for i=tspan(2:end)
    [ti, xi] = ode87(@(t, x) f(t, x, const), [tspan(count),i], x0, options);
    x = [x; xi(2:end,:)]; t = [t; ti(2:end)]; t_idx = [t_idx; numel(t)];
    x0 = x(end,:)'; count = count + 1; 
end

for i=1:length(t)
    x(i,:) = oe2rv(x(i,:), const.mu, 'rad', 'M', 'coe'); 
end

%% Plotting Nominal Trajectories
figure(1); 
plot3(x(:,1),x(:,2), x(:,3), 'r','LineWidth',10,'DisplayName','Nominal');
%scatter3(x(1,1),x(1,2),x(1,3),100,'filled','MarkerFaceColor','r','HandleVisibility','off');
drawnow;
figure(2); 
plot3(x(:,4),x(:,5),x(:,6), 'r','LineWidth', 10, 'DisplayName','Nominal');
%scatter3(x(1,4),x(1,5),x(1,6),100,'filled','MarkerFaceColor','r','HandleVisibility','off');
drawnow;

%% Particle Filter
p.type = 'sharp'; p.color = colors(1,:); p.means = 0; p.display = 0; p.plt = 0; p.name = "PF";  
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/PF/cmake-build-debug/Epochs/Europa/Classical/Revs6";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1; 
for i=[0]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [x_pf, P_pf, n_pf, t_pf(count)] = parse_nongaussian_txt(FILE_PATH, 'coe');

    xest_pf{count} = zeros(size(x_pf(1,:)));
    for j=1:n_pf
        xest_pf{count} = xest_pf{count}+x_pf(j,:).*P_pf(j);
    end

    [~,shp_p] = plot_nongaussian_surface(x_pf(:,1:3),P_pf,[0.14,0.81,0.997],0.5,p);
    [~,shp_v] = plot_nongaussian_surface(x_pf(:,4:6),P_pf,[0.14,0.81,0.997],0.5,p);
    shps_pf{count} = {shp_p, shp_v}; 

    if (i==0)||(i==16)
        p.plt = 1; 
        [axs{axs_count}, ~] = plot_non_gaussian(x_pf,P_pf,t_pf(count),p,axs_count,[0.14,0.81,0.997],0.5,[]);
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end

    count = count + 1;
end

%% UKF (Cartesian)
clear p; p.color = colors(2,:); p.means = 0; p.display = 0; p.plt = 0; p.name = "UKF (Cartesian)";
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/UKF/cmake-build-debug/Epochs/Europa/Cartesian/Revs6";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1; sd = [1,2,3];
for i=[0:16]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt";

    [xest_ukf_c{count}, Sest_ukf_c{count}, t_ukf_c(count)] = parse_gaussian_txt(FILE_PATH, 'cart');
    xesti_ukf = xest_ukf_c{count}; Sesti_ukf = Sest_ukf_c{count};  

    [~,shp_p] = plot_gaussian_ellipsoid(xesti_ukf(1:3,:),Sesti_ukf(1:3,1:3),sd,p);
    [~,shp_v] = plot_gaussian_ellipsoid(xesti_ukf(4:6,:),Sesti_ukf(4:6,4:6),sd,p);
    shps_ukf_c{count} = {shp_p, shp_v}; 

    if (i==0)||(i==16)
        p.plt = 1; 
        plot_gaussian(xest_ukf_c{count},Sest_ukf_c{count},t_ukf_c(count),p,axs_count,sd,axs{axs_count});
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end

    count = count + 1;
end

%% UKF (Equinoctial)
clear p; p.color = colors(3,:); p.means = 0; p.display = 0; p.plt = 0; p.name = "UKF (Equinoctial)";
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/UKF/cmake-build-debug/Epochs/Europa/Equinoctial/Revs6";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1; sd = [1,2,3];
for i=[0]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt";

    [xest_ukf_e{count}, Sest_ukf_e{count}, t_ukf_e(count)] = parse_gaussian_txt(FILE_PATH, 'eoe');
    xesti_ukf = xest_ukf_e{count}; Sesti_ukf = Sest_ukf_e{count};  

    [~,shp_p] = plot_gaussian_ellipsoid(xesti_ukf(1:3,:),Sesti_ukf(1:3,1:3),sd,p);
    [~,shp_v] = plot_gaussian_ellipsoid(xesti_ukf(4:6,:),Sesti_ukf(4:6,4:6),sd,p);
    shps_ukf_e{count} = {shp_p, shp_v}; 

    if (i==0)||(i==16)
        p.plt = 1; 
        plot_gaussian(xest_ukf_e{count},Sest_ukf_e{count},t_ukf_e(count),p,axs_count,sd,axs{axs_count});
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end

    count = count + 1;
end

%% EnKF (Cartesian)
clear p; p.color = colors(4,:); p.means = 0; p.display = 0; p.plt = 0; p.name = "EnKF (Cartesian)";
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Europa/Cartesian/Revs6";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1; sd = [1,2,3];
for i=[0]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt";

    [xest_enkf_c{count}, Sest_enkf_c{count}, t_enkf_c(count)] = parse_gaussian_txt(FILE_PATH, 'cart');
    xesti_enkf = xest_enkf_c{count}; Sesti_enkf = Sest_enkf_c{count};  

    [~,shp_p] = plot_gaussian_ellipsoid(xesti_enkf(1:3,:),Sesti_enkf(1:3,1:3),sd,p);
    [~,shp_v] = plot_gaussian_ellipsoid(xesti_enkf(4:6,:),Sesti_enkf(4:6,4:6),sd,p);
    shps_enkf_c{count} = {shp_p, shp_v}; 

    if (i==0)||(i==16)
        p.plt = 1; 
        plot_gaussian(xest_enkf_c{count},Sest_enkf_c{count},t_enkf_c(count),p,axs_count,sd,axs{axs_count});
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end

    count = count + 1;
end

%% EnKF (Equinoctial)
clear p; p.color = colors(5,:); p.means = 0; p.display = 0; p.plt = 0; p.name = "EnKF (Equinoctial)";
DATA_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Europa/Equinoctial/Revs6";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1; sd = [1,2,3];
for i=[0]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt";

    [xest_enkf_e{count}, Sest_enkf_e{count}, t_enkf_e(count)] = parse_gaussian_txt(FILE_PATH, 'eoe');
    xesti_enkf = xest_enkf_e{count}; Sesti_enkf = Sest_enkf_e{count};  

    [~,shp_p] = plot_gaussian_ellipsoid(xesti_enkf(1:3,:),Sesti_enkf(1:3,1:3),sd,p);
    [~,shp_v] = plot_gaussian_ellipsoid(xesti_enkf(4:6,:),Sesti_enkf(4:6,4:6),sd,p);
    shps_enkf_e{count} = {shp_p, shp_v}; 

    if (i==0)||(i==16)
        p.plt = 1; 
        plot_gaussian(xest_enkf_e{count},Sest_enkf_e{count},t_enkf_e(count),p,axs_count,sd,axs{axs_count});
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end

    count = count + 1;
end

%{
%% GBEES (Equinoctial)
clear p; p.color = colors(6,:); p.plt = 0; p.name = "GBEES (Equinoctial)";
DATA_PATH = "/Users/bhanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/GBEES/GBEES/R2BP/GBEES/cmake-build-debug/Epochs/Europa/Equinoctial/Revs6";
fileList = dir(fullfile(DATA_PATH, '*.txt'));  % List only .txt files
numFiles = numel(fileList);

count = 1; axs_count = 1;
for i=[0]
    FILE_PATH = DATA_PATH + "/pdf_" + num2str(i) + ".txt"; 

    [x_gbees, P_gbees, n_gbees, t_gbees(count)] = parse_nongaussian_txt(FILE_PATH, 'eoe');

    xest_gbees{count} = zeros(size(x_gbees(1,:)));
    for j=1:n_gbees
        xest_gbees{count} = xest_gbees{count}+x_gbees(j,:).*P_gbees(j);
    end

    [~,shp_p] = plot_grid_surface(x_gbees(:,1:3),P_gbees,[0.14,0.81,0.997],0.3,p);
    [~,shp_v] = plot_grid_surface(x_gbees(:,4:6),P_gbees,[0.14,0.81,0.997],0.3,p);
    shps_gbees{count} = {shp_p, shp_v}; 

    if (i==0)||(i==16)
        p.plt = 1; 
        plot_gbees(x_gbees,P_gbees,t_gbees(count),p,axs_count,[0.14,0.81,0.997],0.3,axs{axs_count});
        drawnow; 
        axs_count = axs_count + 1; 
        p.plt = 0;
    end

    count = count + 1;
end
%}

%% Comparison Parameters
% Position/Velocity Error
n = size(xest_pf); n = n(2);
pe = NaN(n,4); ve = NaN(n,4);
for i=1:n
    xi_pf = xest_pf{i}; 
    xi_ukf_c = xest_ukf_c{i};
    xi_ukf_e = xest_ukf_e{i}; 
    xi_enkf_c = xest_enkf_c{i};
    xi_enkf_e = xest_enkf_e{i}; 

    pe(i,1) = norm((x(t_idx(i),1:3)'-xi_ukf_c(1:3)));
    ve(i,1) = norm((x(t_idx(i),4:6)'-xi_ukf_c(4:6)));
    pe(i,2) = norm((x(t_idx(i),1:3)'-xi_ukf_e(1:3)));
    ve(i,2) = norm((x(t_idx(i),4:6)'-xi_ukf_e(4:6)));
    pe(i,3) = norm((x(t_idx(i),1:3)'-xi_enkf_c(1:3)));
    ve(i,3) = norm((x(t_idx(i),4:6)'-xi_enkf_c(4:6)));
    pe(i,4) = norm((x(t_idx(i),1:3)'-xi_enkf_e(1:3)));
    ve(i,4) = norm((x(t_idx(i),4:6)'-xi_enkf_e(4:6)));
end

% Jaccard Coefficient
n = size(xest_pf); n = n(2); 
po_ukf_c = NaN(n,3); vo_ukf_c = NaN(n,3); 
po_ukf_e = NaN(n,3); vo_ukf_e = NaN(n,3); 
po_enkf_c = NaN(n,3); vo_enkf_c = NaN(n,3); 
po_enkf_e = NaN(n,3); vo_enkf_e = NaN(n,3); 
for i=1:n
    shpi_pf = shps_pf{i}; 
    shpi_ukf_c = shps_ukf_c{i}; 
    shpi_ukf_e = shps_ukf_e{i}; 
    shpi_enkf_c = shps_enkf_c{i}; 
    shpi_enkf_e = shps_enkf_e{i}; 

    meani_ukf_c = xest_ukf_c{i}; covi_ukf_c = Sest_ukf_c{i};
    meani_ukf_e = xest_ukf_e{i}; covi_ukf_e = Sest_ukf_e{i};
    meani_enkf_c = xest_enkf_c{i}; covi_enkf_c = Sest_enkf_c{i};
    meani_enkf_e = xest_enkf_e{i}; covi_enkf_e = Sest_enkf_e{i};
    
    pi_pf = shpi_pf{1}; vi_pf = shpi_pf{2}; 
    pi_ukf_c = shpi_ukf_c{1}; vi_ukf_c = shpi_ukf_c{2}; 
    pi_ukf_e = shpi_ukf_e{1}; vi_ukf_e = shpi_ukf_e{2}; 
    pi_enkf_c = shpi_enkf_c{1}; vi_enkf_c = shpi_enkf_c{2}; 
    pi_enkf_e = shpi_enkf_e{1}; vi_enkf_e = shpi_enkf_e{2}; 
    
    for j=1:numel(pi_ukf_c)
        p_ukf_shp.sd = sd(j); 

        p_ukf_shp.Points = pi_ukf_c{j}; 
        p_ukf_shp.mean = meani_ukf_c(1:3); 
        p_ukf_shp.cov = covi_ukf_c(1:3,1:3); 
        po_ukf_c(i,j) = jaccard_coef(pi_pf{j}, p_ukf_shp);
        
        v_ukf_shp.sd = sd(j); 

        v_ukf_shp.Points = vi_ukf_c{j}; 
        v_ukf_shp.mean = meani_ukf_c(4:6); 
        v_ukf_shp.cov = covi_ukf_c(4:6,4:6); 
        vo_ukf_c(i,j) = jaccard_coef(vi_pf{j}, v_ukf_shp);
    end

    for j=1:numel(pi_ukf_e)
        p_ukf_shp.sd = sd(j); 

        p_ukf_shp.Points = pi_ukf_e{j}; 
        p_ukf_shp.mean = meani_ukf_e(1:3); 
        p_ukf_shp.cov = covi_ukf_e(1:3,1:3); 
        po_ukf_e(i,j) = jaccard_coef(pi_pf{j}, p_ukf_shp);
        
        v_ukf_shp.sd = sd(j); 

        v_ukf_shp.Points = vi_ukf_e{j}; 
        v_ukf_shp.mean = meani_ukf_e(4:6); 
        v_ukf_shp.cov = covi_ukf_e(4:6,4:6); 
        vo_ukf_e(i,j) = jaccard_coef(vi_pf{j}, v_ukf_shp);
    end
    
    for j=1:numel(pi_enkf_c)
        p_enkf_shp.sd = sd(j); 

        p_enkf_shp.Points = pi_enkf_c{j}; 
        p_enkf_shp.mean = meani_enkf_c(1:3); 
        p_enkf_shp.cov = covi_enkf_c(1:3,1:3); 
        po_enkf_c(i,j) = jaccard_coef(pi_pf{j}, p_enkf_shp);
        
        v_enkf_shp.sd = sd(j); 

        v_enkf_shp.Points = vi_enkf_c{j}; 
        v_enkf_shp.mean = meani_enkf_c(4:6); 
        v_enkf_shp.cov = covi_enkf_c(4:6,4:6); 
        vo_enkf_c(i,j) = jaccard_coef(vi_pf{j}, v_enkf_shp);
    end

    for j=1:numel(pi_enkf_e)
        p_enkf_shp.sd = sd(j); 

        p_enkf_shp.Points = pi_enkf_e{j}; 
        p_enkf_shp.mean = meani_enkf_e(1:3); 
        p_enkf_shp.cov = covi_enkf_e(1:3,1:3); 
        po_enkf_e(i,j) = jaccard_coef(pi_pf{j}, p_enkf_shp);
        
        v_enkf_shp.sd = sd(j); 

        v_enkf_shp.Points = vi_enkf_e{j}; 
        v_enkf_shp.mean = meani_enkf_e(4:6); 
        v_enkf_shp.cov = covi_enkf_e(4:6,4:6); 
        vo_enkf_e(i,j) = jaccard_coef(vi_pf{j}, v_enkf_shp);
    end
end

% Timing
FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/PF/cmake-build-debug/Epochs/Europa/Classical/timing6.txt";
[pt_pf, st_pf] = parse_timing_txt(FILE_PATH);

FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/UKF/cmake-build-debug/Epochs/Europa/Equinoctial/timing6.txt";
[pt_ukf_e, st_ukf_e] = parse_timing_txt(FILE_PATH);
FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/UKF/cmake-build-debug/Epochs/Europa/Cartesian/timing6.txt";
[pt_ukf_c, st_ukf_c] = parse_timing_txt(FILE_PATH);

FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Europa/Equinoctial/timing6.txt";
[pt_enkf_e, st_enkf_e] = parse_timing_txt(FILE_PATH);
FILE_PATH = "/Users/bhanson/OneDrive - UC San Diego/UCSD/Research/GBEES/GBEES/R2BP/EnKF/cmake-build-debug/Epochs/Europa/Cartesian/timing6.txt";
[pt_enkf_c, st_enkf_c] = parse_timing_txt(FILE_PATH);

pt_ukf_e = pt_ukf_e./pt_pf; pt_ukf_e(1) = 1; 
pt_ukf_c = pt_ukf_c./pt_pf; pt_ukf_c(1) = 1; 
pt_enkf_e = pt_enkf_e./pt_pf; pt_enkf_e(1) = 1; 
pt_enkf_c = pt_enkf_c./pt_pf; pt_enkf_c(1) = 1; 

fnew = figure; hold on; fnew.Position = [50 150 1400 800]; 
tiledlayout(3, 1, "TileSpacing",'compact');

ax1 = nexttile; hold on;legend('Location','northwest','Interpreter','latex');
plot(t_pf,pe(:,1),'-o', 'Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Cartesian)')
plot(t_pf,pe(:,2),'-o', 'Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Equinoctial)')
plot(t_pf,pe(:,3),'-o', 'Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Cartesian)')
plot(t_pf,pe(:,4),'-o', 'Color', colors(5,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Equinoctial)')

ax2 = nexttile; hold on;legend('Location','southwest','Interpreter','latex');
plot(t_pf,po_ukf_c(:,1),'-o','Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Cartesian) - 1$\sigma$')
plot(t_pf,po_ukf_c(:,2),'-square','Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Cartesian) - 2$\sigma$')
plot(t_pf,po_ukf_c(:,3),'-diamond','Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Cartesian) - 3$\sigma$')
plot(t_pf,po_ukf_e(:,1),'-o','Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Equinoctial) - 1$\sigma$')
plot(t_pf,po_ukf_e(:,2),'-square','Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Equinoctial) - 2$\sigma$')
plot(t_pf,po_ukf_e(:,3),'-diamond','Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Equinoctial) - 3$\sigma$')
plot(t_pf,po_enkf_c(:,1),'-o','Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Cartesian) - 1$\sigma$')
plot(t_pf,po_enkf_c(:,2),'-square','Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Cartesian) - 2$\sigma$')
plot(t_pf,po_enkf_c(:,3),'-diamond','Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Cartesian) - 3$\sigma$')
plot(t_pf,po_enkf_e(:,1),'-o','Color', colors(5,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Equinoctial) - 1$\sigma$')
plot(t_pf,po_enkf_e(:,2),'-square','Color', colors(5,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Equinoctial) - 2$\sigma$')
plot(t_pf,po_enkf_e(:,3),'-diamond','Color', colors(5,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Equinoctial) - 3$\sigma$')

ax3 = nexttile; hold on;legend('Location','northwest','Interpreter','latex');
plot(st_ukf_c,pt_ukf_c,'-o', 'Color', colors(2,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Cartesian)')
plot(st_ukf_e,pt_ukf_e,'-o', 'Color', colors(3,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','UKF (Equinoctial)')
plot(st_enkf_c,pt_enkf_c,'-o', 'Color', colors(4,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Cartesian)')
plot(st_enkf_e,pt_enkf_e,'-o', 'Color', colors(5,:), 'LineWidth', 1, 'MarkerSize', 10, 'DisplayName','EnKF (Equinoctial)')

tickvals = [0:0.25:round(t_pf(end)/const.T)];
ticklabels = {}; count = 1;  for i=tickvals, ticklabels{count} = num2str(i); count = count + 1; end

set(ax1 , ...
     'FontName' , 'Times' , ...
     'FontSize' , 10 , ... 
     'XLim', [t_pf(1),1.001*t_pf(end)], ...
     'TickDir' , 'out' , ...
     'TickLength' , [0.0125 0.0125] , ...
     'XTick', [0:0.25*const.T:round(t_pf(end)/const.T)*const.T], ...
     'XTickLabel', ticklabels, ...
     'XGrid' , 'on' , 'YGrid' , 'on' )  
ylabel( ax1 , "Position Error (km)", ...
                'Interpreter', 'LaTex' , 'Rotation', 90 ); 
ax1.YLabel.HorizontalAlignment = 'center';          
set(ax2 , ...
     'FontName' , 'Times' , ...
     'FontSize' , 10 , ... 
     'XLim', [t_pf(1),1.001*t_pf(end)], ...
     'TickDir' , 'out' , ...
     'TickLength' , [0.0125 0.0125] , ...
     'XTick', [0:0.25*const.T:round(t_pf(end)/const.T)*const.T], ...
     'XTickLabel', ticklabels, ...
     'XGrid' , 'on' , 'YGrid' , 'on' ) 
ylabel(ax2, "$J_{\mathbf{r}}$", ...
                'Interpreter', 'LaTex' , 'Rotation', 90 ); 
ax2.YLabel.HorizontalAlignment = 'center';    
set(ax3 , ...
     'FontName' , 'Times' , ...
     'FontSize' , 10 , ... 
     'XLim', [t_pf(1),1.001*t_pf(end)], ...
     'YLim', [0,1.1], ...
     'TickDir' , 'out' , ...
     'TickLength' , [0.0125 0.0125] , ...
     'XTick', [0:0.25*const.T:round(t_pf(end)/const.T)*const.T], ...
     'XTickLabel', ticklabels, ...
     'XGrid' , 'on' , 'YGrid' , 'on' ) 
xlabel(ax3, "$t$ (orbital periods)", ...
                'Interpreter', 'LaTex' , 'Rotation', 0 ); 
ylabel(ax3, "Normalized computation time", ...
                'Interpreter', 'LaTex' , 'Rotation', 90 ); 
ax3.YLabel.HorizontalAlignment = 'center';   
%}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = f(t, x,const)
    x1 = [0; 0; 0; 0; 0; sqrt(const.mu/x(1)^3)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axs,shps] = plot_non_gaussian(x,P,t,p,count,isovalues,alpha,axs)
    if isempty(axs)
        p.display = 0; 

        figure(1); 
        [~,shp_p] = plot_nongaussian_surface(x(:,1:3),P,isovalues,alpha,p);
        drawnow;
        figure(2);
        [~,shp_v] = plot_nongaussian_surface(x(:,4:6),P,isovalues,alpha,p);
        drawnow;

        shps = {shp_p, shp_v}; 
    
        fi = figure(2*count+1); hold on; fi.Position = [50 100 700 700];  
        set(gca, 'FontName' , 'Times','FontSize',12); 
        tiledlayout(2,2);
    
        nexttile;  hold on;
        plot_nongaussian_surface(x(:,1:2),P,isovalues,alpha,p);
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile; hold on; 
        plot_nongaussian_surface(x(:,[3,2]),P,isovalues,alpha,p);
        xlabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on;
        plot_nongaussian_surface(x(:,[1,3]),P,isovalues,alpha,p);
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; lighting phong; light('Position',[1 -1 1]);
        p.display=1; 
        plot_nongaussian_surface(x(:,1:3),P,isovalues,alpha,p); 
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        zlabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        view(45, 30); 
    
        hL = legend('Orientation', 'horizontal');
        hL.Layout.Tile = 'north';
        sgtitle(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize',18, 'FontName', 'Times', 'Interpreter','latex');
        axs_p = {nexttile(1),nexttile(2),nexttile(3),nexttile(4)};
    
        p.display=0; 
        fj = figure(2*count+2); hold on; fj.Position = [800 100 700 700];
        set(gca, 'FontName' , 'Times','FontSize',12); 
        tiledlayout(2,2);
    
        nexttile;  hold on;
        plot_nongaussian_surface(x(:,[4,5]),P,isovalues,alpha,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; 
        plot_nongaussian_surface(x(:,[6,5]),P,isovalues,alpha,p);
        xlabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on;
        plot_nongaussian_surface(x(:,[4,6]),P,isovalues,alpha,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; lighting phong; light('Position',[1 -1 1]);
        p.display=1; 
        plot_nongaussian_surface(x(:,4:6),P,isovalues,alpha,p);  
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        zlabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        view(45, 30); 
    
        hL = legend('Orientation', 'horizontal');
        hL.Layout.Tile = 'north';
        sgtitle(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize',18, 'FontName', 'Times', 'Interpreter','latex');
        axs_v={nexttile(1),nexttile(2),nexttile(3),nexttile(4)};

        axs = {axs_p,axs_v};
    else
        p.display = 0; axs_p = axs{1}; axs_v = axs{2}; 

        figure(1);
        [~,shp_p] = plot_nongaussian_surface(x(:,1:3),P,isovalues,alpha,p);
        drawnow;
        figure(2);
        [~,shp_v] = plot_nongaussian_surface(x(:,4:6),P,isovalues,alpha,p);
        drawnow;
        
        shps = {shp_p,shp_v};

        p.axh = axs_p{1}; 
        plot_nongaussian_surface(x(:,1:2),P,isovalues,alpha,p);

        p.axh = axs_p{2}; 
        plot_nongaussian_surface(x(:,[3,2]),P,isovalues,alpha,p);

        p.axh = axs_p{3}; 
        plot_nongaussian_surface(x(:,[1,3]),P,isovalues,alpha,p);
    
        p.axh = axs_p{4}; p.display=1;
        plot_nongaussian_surface(x(:,1:3),P,isovalues,alpha,p);
        p.display=0;
        
        p.axh = axs_v{1}; 
        plot_nongaussian_surface(x(:,4:5),P,isovalues,alpha,p);
        
        p.axh = axs_v{2}; 
        plot_nongaussian_surface(x(:,[6,5]),P,isovalues,alpha,p);
        
        p.axh = axs_v{3}; 
        plot_nongaussian_surface(x(:,[4,6]),P,isovalues,alpha,p);
    
        p.axh = axs_v{4}; p.display=1;
        plot_nongaussian_surface(x(:,4:6),P,isovalues,alpha,p);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axs,shps] = plot_gbees(x,P,t,p,count,isovalues,alpha,axs)
    if isempty(axs)
        p.display = 0; 

        figure(1); 
        [~,shp_p] = plot_grid_surface(x(:,1:3),P,isovalues,alpha,p);
        drawnow;
        figure(2);
        [~,shp_v] = plot_grid_surface(x(:,4:6),P,isovalues,alpha,p);
        drawnow;

        shps = {shp_p, shp_v}; 
    
        fi = figure(2*count+1); hold on; fi.Position = [50 100 700 700];  
        set(gca, 'FontName' , 'Times','FontSize',12); 
        tiledlayout(2,2);
    
        nexttile;  hold on;
        plot_grid_surface(x(:,1:2),P,isovalues,alpha,p);
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile; hold on; 
        plot_grid_surface(x(:,[3,2]),P,isovalues,alpha,p);
        xlabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on;
        plot_grid_surface(x(:,[1,3]),P,isovalues,alpha,p);
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; lighting phong; light('Position',[1 -1 1]);
        p.display=1; 
        plot_grid_surface(x(:,1:3),P,isovalues,alpha,p); 
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        zlabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        view(45, 30); 
    
        hL = legend('Orientation', 'horizontal');
        hL.Layout.Tile = 'north';
        sgtitle(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize',18, 'FontName', 'Times', 'Interpreter','latex');
        axs_p = {nexttile(1),nexttile(2),nexttile(3),nexttile(4)};
    
        p.display=0; 
        fj = figure(2*count+2); hold on; fj.Position = [800 100 700 700];
        set(gca, 'FontName' , 'Times','FontSize',12); 
        tiledlayout(2,2);
    
        nexttile;  hold on;
        plot_grid_surface(x(:,[4,5]),P,isovalues,alpha,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; 
        plot_grid_surface(x(:,[6,5]),P,isovalues,alpha,p);
        xlabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on;
        plot_grid_surface(x(:,[4,6]),P,isovalues,alpha,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; lighting phong; light('Position',[1 -1 1]);
        p.display=1; 
        plot_grid_surface(x(:,4:6),P,isovalues,alpha,p);  
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        zlabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        view(45, 30); 
    
        hL = legend('Orientation', 'horizontal');
        hL.Layout.Tile = 'north';
        sgtitle(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize',18, 'FontName', 'Times', 'Interpreter','latex');
        axs_v={nexttile(1),nexttile(2),nexttile(3),nexttile(4)};

        axs = {axs_p,axs_v};
    else
        p.display = 0; axs_p = axs{1}; axs_v = axs{2}; 

        figure(1);
        [~,shp_p] = plot_grid_surface(x(:,1:3),P,isovalues,alpha,p);
        drawnow;
        figure(2);
        [~,shp_v] = plot_grid_surface(x(:,4:6),P,isovalues,alpha,p);
        drawnow;
        
        shps = {shp_p,shp_v};

        p.axh = axs_p{1}; 
        plot_grid_surface(x(:,1:2),P,isovalues,alpha,p);

        p.axh = axs_p{2}; 
        plot_grid_surface(x(:,[3,2]),P,isovalues,alpha,p);

        p.axh = axs_p{3}; 
        plot_grid_surface(x(:,[1,3]),P,isovalues,alpha,p);
    
        p.axh = axs_p{4}; p.display=1;
        plot_grid_surface(x(:,1:3),P,isovalues,alpha,p);
        p.display=0;
        
        p.axh = axs_v{1}; 
        plot_grid_surface(x(:,4:5),P,isovalues,alpha,p);
        
        p.axh = axs_v{2}; 
        plot_grid_surface(x(:,[6,5]),P,isovalues,alpha,p);
        
        p.axh = axs_v{3}; 
        plot_grid_surface(x(:,[4,6]),P,isovalues,alpha,p);
    
        p.axh = axs_v{4}; p.display=1;
        plot_grid_surface(x(:,4:6),P,isovalues,alpha,p);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [axs, shps] = plot_gaussian(x,S,t,p,count,sd,axs)
    if isempty(axs)
        p.display = 0;
    
        figure(1);
        [~, shp_p] = plot_gaussian_ellipsoid(x(1:3),S(1:3,1:3),sd,p);
        drawnow;
        figure(2);
        [~, shp_v] = plot_gaussian_ellipsoid(x(4:6),S(4:6,4:6),sd,p);
        drawnow;

        shps = {shp_p,shp_v};

        fi = figure(2*count+1); hold on; fi.Position = [50 100 700 700];  
        set(gca, 'FontName' , 'Times','FontSize',12); 
        tiledlayout(2,2);
    
        nexttile;  hold on;
        plot_gaussian_ellipsoid(x(1:2),S(1:2,1:2),sd,p);   
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile; hold on; 
        pplot_gaussian_ellipsoid(x([3,2]),S([3,2],[3,2]),sd,p);
        xlabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on;
        plot_gaussian_ellipsoid(x([1,3]),S([1,3],[1,3]),sd,p);
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; lighting phong; light('Position',[1 -1 1]);
        p.display=1; 
        plot_gaussian_ellipsoid(x(1:3),S(1:3,1:3),sd,p);
        xlabel('x (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('y (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        zlabel('z (km)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        view(45, 30); 
    
        hL = legend('Orientation', 'horizontal');
        hL.Layout.Tile = 'north';
        sgtitle(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize',18, 'FontName', 'Times', 'Interpreter','latex');
        axs_p = {nexttile(1),nexttile(2),nexttile(3),nexttile(4)};
    
        p.display=0; 
        fj = figure(2*count+2); hold on; fj.Position = [800 100 700 700];
        set(gca, 'FontName' , 'Times','FontSize',12); 
        tiledlayout(2,2);
    
        nexttile;  hold on;
        plot_gaussian_ellipsoid(x(4:5),S(4:5,4:5),sd,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; 
        plot_gaussian_ellipsoid(x([6,5]),S([6,5],[6,5]),sd,p);
        xlabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on;
        plot_gaussian_ellipsoid(x([4,6]),S([4,6],[4,6]),sd,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
    
        nexttile;  hold on; lighting phong; light('Position',[1 -1 1]);
        p.display=1; 
        plot_gaussian_ellipsoid(x(4:6),S(4:6,4:6),sd,p);
        xlabel('$v_x$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        ylabel('$v_y$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        zlabel('$v_z$ (km/s)', 'FontSize',12, 'FontName', 'Times', 'Interpreter','latex');
        view(45, 30); 
    
        hL = legend('Orientation', 'horizontal');
        hL.Layout.Tile = 'north';
        sgtitle(append('Orbit time: ', num2str(t/3600), ' hrs'), 'FontSize',18, 'FontName', 'Times', 'Interpreter','latex');
        axs_v={nexttile(1),nexttile(2),nexttile(3),nexttile(4)};

        axs = {axs_p,axs_v};
    else
        p.display = 0; axs_p = axs{1}; axs_v = axs{2};  
    
        figure(1);
        [~, shp_p] = plot_gaussian_ellipsoid(x(1:3),S(1:3,1:3),sd,p);
        drawnow;
        figure(2);
        [~, shp_v] = plot_gaussian_ellipsoid(x(4:6),S(4:6,4:6),sd,p);
        drawnow;

        shps = {shp_p,shp_v};

        p.axh = axs_p{1}; 
        plot_gaussian_ellipsoid(x(1:2),S(1:2,1:2),sd,p);    
        
        p.axh = axs_p{2}; 
        plot_gaussian_ellipsoid(x([3,2]),S([3,2],[3,2]),sd,p);
        
        p.axh = axs_p{3}; 
        plot_gaussian_ellipsoid(x([1,3]),S([1,3],[1,3]),sd,p);
    
        p.axh = axs_p{4}; p.display=1;
        plot_gaussian_ellipsoid(x(1:3),S(1:3,1:3),sd,p);
        p.display=0;

        p.axh = axs_v{1}; 
        plot_gaussian_ellipsoid(x(4:5),S(4:5,4:5),sd,p);
        
        p.axh = axs_v{2}; 
        plot_gaussian_ellipsoid(x([6,5]),S([6,5],[6,5]),sd,p);

        p.axh = axs_v{3}; 
        plot_gaussian_ellipsoid(x([4,6]),S([4,6],[4,6]),sd,p);
    
        p.axh = axs_v{4}; p.display=1;
        plot_gaussian_ellipsoid(x(4:6),S(4:6,4:6),sd,p);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, P, n, t] = parse_nongaussian_txt(filename,oetype)
    fileID = fopen(filename, 'r'); t = str2double(fgetl(fileID));
    
    count = 1; 
    while ~feof(fileID)
        line = split(fgetl(fileID)); % Read a line as a string
        P(count,1) = str2double(line{1});
        oe = [str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5});str2double(line{6});str2double(line{7})];
        if strcmp(oetype, 'cart')
            x(count, :) = oe;
        elseif strcmp(oetype, 'coe')    
            x(count, :) = oe2rv(oe,3202.738774922892,'rad','M',oetype); % Convert to Cartesian
        elseif strcmp(oetype, 'eoe')
            x(count, :) = oe2rv(oe,3202.738774922892,'rad',[],oetype); % Convert to Cartesian
        end
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
    
    if strcmp(oetype, 'cart')
        x = oe;
    elseif strcmp(oetype, 'coe')    
        x = oe2rv(oe,3202.738774922892,'rad','M',oetype); % Convert to Cartesian
    elseif strcmp(oetype, 'eoe')
        x = oe2rv(oe,3202.738774922892,'rad',[],oetype); % Convert to Cartesian
    end

    fgetl(fileID); % skip blank line

    for i=1:6
        line = split(fgetl(fileID)); % Read a line as a string
        S(i,:) = [str2double(line{1});str2double(line{2});str2double(line{3});str2double(line{4});str2double(line{5});str2double(line{6})];
    end
    
    if strcmp(oetype, 'coe')  
        S = oe2rv_cov(oe,S,3202.738774922892,'M',oetype); % Convert to Cartesian
    elseif strcmp(oetype, 'eoe')
        S = oe2rv_cov(oe,S,3202.738774922892,[],oetype); % Convert to Cartesian 
    end

    % Close the file
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pt, st] = parse_timing_txt(FILE_PATH)
    fileID = fopen(FILE_PATH, 'r');
    pt = []; st = [];
    
    while ~feof(fileID)
        line = split(fgetl(fileID)); 
        pt(end+1) = str2double(line{3});
        st(end+1) = str2double(line{7});
    end
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%