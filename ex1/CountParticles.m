function C=CountParticles(pos,L,M)
% Author: Giorgio Palermo
% V3 .0
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
% to test:
% 
% clear
% 
% N = 30; %number of particles
% L=1; %size of the box
% M = 10; %number of voxels
% pos = L*rand(N,3);
% 
% CountParticles(pos,L,M);

dx = L/M;
C=zeros(M,M,M);

for i=1:size(pos,1)
    idx = round(pos(i,1)/dx + .5);
    idy = round(pos(i,2)/dx + .5);
    idz = round(pos(i,3)/dx + .5);
    C(idx, idy, idz) = C(idx,idy,idz) + 1;
end
end