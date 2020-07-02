% 2D Random Walk plot

dt = 1e-5; %timestep, seconds
T = 200*dt; %s, duration of the simulation
% N = round(T/dt + .5);  % number of steps
D = 200;    %um/s, diffusion coefficient

x=diffuse3(D,dt,T);
figure(1)
scatter3(x(:,1),x(:,2),x(:,3));