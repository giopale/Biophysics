% CountParticles 
% Author: Giorgio Palermo
% V1.0
% 
% 
% Syntax C=CountParticles(pos,L,M);
% 
% Input:
% pos is a Nx3 array containing the positions of N
% particles randomly distributed within the box.
% L is the box side.
% 
% M is the number of voxels along one dimension.
% 
% Output:
% C is a MxMxM matrix having the (i,j,k)-th element
% equal to the number of particles within the (i,j,k)- 
% th voxel of the box.
clear

N = 30; %number of particles
L=1; %size of the box
M = 10; %number of voxels
pos = L*rand(N,1);

Y=zeros(N,1);
Y(:,1) = discretize(pos(:,1),M);
Z=Y - ceil(pos/M)

