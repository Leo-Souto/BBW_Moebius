function [V,F] = adaptive_mesh(P_handles)

TR = IcosahedronMesh;
TR = SubdivideSphericalMesh(TR, 4);
P_ico = TR.Points;
size(P_ico)
distances = zeros(size(P_ico,1),size(P_handles,1));
for i=1:size(P_ico,1)
    for j = 1: size(P_handles,1)
        distances(i,j) = norm(P_ico(i,:)-P_handles(j,:));
    end
end
distances = min(distances,[],2);
probabilities = (1-distances/(max(distances)));

V = P_ico(rand(size(distances))<probabilities.^2,:);
F = convhull(V);
size(V)

tsurf(F,V); axis equal; cameratoolbar;

% [VV,~,FF] = tetgen(V,F,'Flags', '-q1.0/8');
% figure
% tsurf(FF,VV); axis equal; cameratoolbar;
% 
% 
% V = VV;
% F = FF;