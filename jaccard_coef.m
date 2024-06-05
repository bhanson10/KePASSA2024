function J = jaccard_coef(shp1, shp2, npts) 
% jaccard_coef.m
% Benjamin Hanson, 2024
% 
% Given two 3D isosurface shapes (either alphaShape for a non-Gaussian
% distribution or ellipse for a Gaussian distribution), approximate the
% Jaccard coefficient 
% 
% Inputs:
%          shp1 -- 2D/3D shape 1
%          shp2 -- 2D/3D shape 2
%          npts -- Number of points in Monte Carlo approximation
%
% Outputs:
%          overlap -- Volume overlap of the two shapes

if ~exist('npts', 'var'), npts = 50; end
N = size(shp1.Points); 

switch N(2)
    case 2, J = jaccard_coef2D(shp1, shp2, npts); 
    case 3, J = jaccard_coef3D(shp1, shp2, npts); 
   otherwise
      error('Unsupported dimensionality');
end
end

function J = jaccard_coef2D(shp1, shp2, npts) 
isAlpha1 = isa(shp1, 'alphaShape'); 
isAlpha2 = isa(shp2, 'alphaShape'); 

x_min = min(shp1.Points(:,1)); x_max = max(shp1.Points(:,1));
y_min = min(shp1.Points(:,2)); y_max = max(shp1.Points(:,2));

x_min = min(x_min,min(shp2.Points(:,1))); x_max = max(x_max,max(shp2.Points(:,1)));
y_min = min(y_min,min(shp2.Points(:,2))); y_max = max(y_max,max(shp2.Points(:,2)));

count_union = 0; count_intersection = 0; 
for i=linspace(x_min,x_max,npts)
    for j=linspace(y_min,y_max,npts)
        pt = [i, j]; 

        % Checking in Shp1
        inShp1 = 0; 
        if isAlpha1
            if inShape(shp1,pt) 
                inShp1 = 1; 
            end
        else
            dM = sqrt((pt' - shp1.mean)'*inv(shp1.cov)*(pt' - shp1.mean));
            if dM <= shp1.sd
                inShp1 = 1; 
            end
        end

        % Checking in Shp2
        inShp2 = 0; 
        if isAlpha2
            if inShape(shp2,pt) 
                inShp2 = 1; 
            end
        else
            dM = sqrt((pt' - shp2.mean)'*inv(shp2.cov)*(pt' - shp2.mean));
            if dM <= shp2.sd
                inShp2 = 1; 
            end
        end

        if inShp1||inShp2
            count_union = count_union + 1; 
        end

        if inShp1&&inShp2
            count_intersection = count_intersection + 1; 
        end
    end
end

J = count_intersection / count_union; 
end

function J = jaccard_coef3D(shp1, shp2, npts) 
isAlpha1 = isa(shp1, 'alphaShape'); 
isAlpha2 = isa(shp2, 'alphaShape'); 

x_min = min(shp1.Points(:,1)); x_max = max(shp1.Points(:,1));
y_min = min(shp1.Points(:,2)); y_max = max(shp1.Points(:,2));
z_min = min(shp1.Points(:,3)); z_max = max(shp1.Points(:,3));

x_min = min(x_min,min(shp2.Points(:,1))); x_max = max(x_max,max(shp2.Points(:,1)));
y_min = min(y_min,min(shp2.Points(:,2))); y_max = max(y_max,max(shp2.Points(:,2)));
z_min = min(z_min,min(shp2.Points(:,3))); z_max = max(z_max,max(shp2.Points(:,3)));

count_union = 0; count_intersection = 0; 
for i=linspace(x_min,x_max,npts)
    for j=linspace(y_min,y_max,npts)
        for k=linspace(z_min,z_max,npts)
            pt = [i, j, k]; 

            % Checking in Shp1
            inShp1 = 0; 
            if isAlpha1
                if inShape(shp1,pt) 
                    inShp1 = 1; 
                end
            else
                dM = sqrt((pt' - shp1.mean)'*inv(shp1.cov)*(pt' - shp1.mean));
                if dM <= shp1.sd
                    inShp1 = 1; 
                end
            end

            % Checking in Shp2
            inShp2 = 0; 
            if isAlpha2
                if inShape(shp2,pt) 
                    inShp2 = 1; 
                end
            else
                dM = sqrt((pt' - shp2.mean)'*inv(shp2.cov)*(pt' - shp2.mean));
                if dM <= shp2.sd
                    inShp2 = 1; 
                end
            end

            if inShp1||inShp2
                count_union = count_union + 1; 
            end

            if inShp1&&inShp2
                count_intersection = count_intersection + 1; 
            end
        end
    end
end

J = count_intersection / count_union; 
end