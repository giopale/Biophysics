function [C,Npos]=CountParticles2(pos,L,M)
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
% CountParticles(pos,L,M);

dx = L/M;
C=zeros(M,M,M);
N=size(pos,1);

idx= zeros(N,3);
for i=1:N
    
    idx(i,1) = round(pos(i,1)/dx + .5);
    idx(i,2) = round(pos(i,2)/dx + .5);
    idx(i,3) = round(pos(i,3)/dx + .5);
    C(idx(i,1), idx(i,2), idx(i,3)) = C(idx(i,1), idx(i,2), idx(i,3)) + 1;
    
end
for i=1:N
    Npos(i,1) = C(idx(i,1), idx(i,2), idx(i,3));
end
% [~,~,ic] = unique(idx,'rows','stable');
% h=accumarray(ic,1);
% Npos = h(ic);

% T=table(pos,idx,  ic, Npos)
% h
end










