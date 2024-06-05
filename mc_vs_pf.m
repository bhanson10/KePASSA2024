clear all; clc; close all; 

initialize_figures('n', 1:2, 'bg', 'clear', 'spacing', {[251,84,903,782]});

mu = [0; 0]; P = [1 0; 0 1]; 

n = 100000;

x = mvnrnd(mu,P,n);
P = gd(x,mu,P,n); P = P./sum(P);

figure(1) 
scatter(x(:,1), x(:,2), 10, 'k', 'filled')

figure(2)
colors = jet(100);
scaled_P = (P - min(P))./(max(P)-min(P));
colors = colors(round(scaled_P.*99)+1,:);
scatter(x(:,1), x(:,2), 10, colors, 'filled')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = gd(x,mu,P,n)
    p = zeros(n,1); 

    for i=1:n
        p(i) = exp(-0.5*((x(i,:)'-mu)'*inv(P)*(x(i,:)'-mu))); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
