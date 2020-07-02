function [C,Npos]=CountParticlesBest(pos,L,M)
% Author: Giorgio Palermo
% V4 .0
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
%
% Npos is a Nx1 array having the i-th element equal 
% to the number of particles within the voxel of the 
% i-th particle in pos matrix.
% 
% to test:
% 
% clear
% 
% N = 30; %number of particles
% L=1; %size of the box
% M = 10; %number of voxels
% pos = L*rand(N,3);
% 
% CountParticlesBest(pos,L,M);

C=zeros(M,M,M,'uint32');
N=size(pos,1);
Npos = zeros(N,1,'uint32');
scale_factor = M/L;

pos = ceil(pos*scale_factor);

loc = sub2ind([M M M],pos(:,1),pos(:,2),pos(:,3));
for k=1:N
    C(loc(k)) = C(loc(k))+1;
end

Npos = C(loc);










