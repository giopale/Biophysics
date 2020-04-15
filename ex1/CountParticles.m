% CountParticles 
% Author: Giorgio Palermo
% V2.0
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

function C=CountParticles(pos,M)
Y(:,1) = discretize(pos(:,1),M);
Y(:,2) = discretize(pos(:,2),M);
Y(:,3) = discretize(pos(:,3),M);
%disp(Y);
C = zeros(M,M,M);
 for i=1:size(Y)
     idx = sub2ind(size(C),Y(i,1),Y(i,2),Y(i,3));
     C(idx)= C(idx) + 1;
 end
%     disp(C) 
    
end