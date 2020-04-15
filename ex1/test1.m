
clear

N = 1e6; %number of particles
L=1; %size of the box
M = 100; %number of voxels
pos = L*randn(N,3);

tic
C =CountParticles(pos,M);
toc

