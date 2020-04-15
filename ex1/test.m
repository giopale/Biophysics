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

N = 1e6; %number of particles
L=1; %size of the box
M = 100; %number of voxels
pos = L*randn(N,3);

tic
C =CountParticles(pos,L,M);
toc

function C=CountParticles(pos,L,M)
C = zeros(M,M,M);

    for i=1:M
    C(:,:,i) = histcounts2(pos(:,1),pos(:,2),M);
    end
end