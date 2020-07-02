
clear

if 1<2
    N = 1e7; %number of particles
    L=10; %size of the box
    M = 300; %number of voxels
    pos = L*rand(N,3);

else
    N = 30; %number of particles
    L=10; %size of the box
    M = 2; %number of voxels
    pos = L*rand(N,3);
end

tic
[C,Npos] = CountParticlesBest(pos,L,M);
toc

