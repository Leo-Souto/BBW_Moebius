function [V,F] = adaptive_mesh(P_handles,max_iter)
% ADAPTIVE_MESH Computes a mesh that is finer close to handles
%
% [V,F] = adaptive_mesh(P_handles,max_iter)
%
% Inputs:
%   P_handles: #P_handles x 3 matrix of handle coordinates
%   max_iter: maximum number of subdivisions
% Output:
%   V: #V x 3 matrix of 3d vertices
%   F: #F x 3 matrix of indices into V

% initial mesh
TR = IcosahedronMesh;
TR = SubdivideSphericalMesh(TR, 1);

for iter=1:max_iter

    % barycenters of each face
    B = barycenter(TR.Points,TR.ConnectivityList);

    % for each barycenter, determine the smallest (geodesic) distance to a
    % handle point
    distances = acos(B*P_handles');
    distances = min(distances,[],2);

    % apply kumaraswamy function to distances
    kumar = kumaraswamy(distances,1,max_iter,pi,0);
    
    % choose the triangles with higher kumaraswamy function values to
    % subdivide once more
    NewConnectivityList = TR.ConnectivityList(kumar>=iter,:);
    if (~isempty(NewConnectivityList))
        TR_new =  triangulation(NewConnectivityList,TR.Points);
        TR_sub = SubdivideSphericalMesh(TR_new, 1);
    end

    % append new points to proveious point list and remove duplicates
    if (~isempty(NewConnectivityList))
        final_points = [TR.Points; TR_sub.Points];
    else
        final_points = [TR.Points];
    end
    [final_points,~,~] = remove_duplicate_vertices(final_points,1e-7);
    
    % calculate the convex hull of the final points
    TR_conv = convhulln(final_points);

    % update TR
    clear TR;
    TR.Points = final_points;
    TR.ConnectivityList = TR_conv;
end

% return values
F = TR_conv;
V = final_points;
