function [ sphere voxels ] = mvpaa_makeSphere( radius )
%MAKE_SPHERE Makes a sphere
% The sphere is centred on one central voxel.
% The radius excludes the central voxel, so diameter is radius*2 + 1.

sphere=zeros(ceil(radius)*2 + 1, ceil(radius)*2 + 1,ceil(radius)*2 + 1);
[X Y Z] = meshgrid(1:(ceil(radius)*2 + 1), 1:(ceil(radius)*2 + 1), 1:(ceil(radius)*2 + 1));
X = X - (ceil(radius) + 1);
Y = Y - (ceil(radius) + 1);
Z = Z - (ceil(radius) + 1);
D = sqrt(X.^2 + Y.^2 + Z.^2);

sphere(D < (radius + 0.001)) = 1;
voxels = sum(sum(sum(sphere)));