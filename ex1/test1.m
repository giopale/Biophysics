
clear

if 1<2 
    N = 1e6; %number of particles
    L=10; %size of the box
    M = 100; %number of voxels
    pos = L*rand(N,3);

else
    N = 30; %number of particles
    L=10; %size of the box
    M = 2; %number of voxels
    pos = L*rand(N,3);
end

tic
[C,Npos] = CountParticles2(pos,L,M);
toc

