function [V,F] = adaptive_mesh(P_handles)

% NEW CODE (CURRENTLY REFINES TRIANGLES WITH ALL VERTICES WITH Z-COORDINATE
% GRATER THAN ZERO)
TR = IcosahedronMesh;
TR = SubdivideSphericalMesh(TR, 1);
faces = [];
for k=1:size(TR.ConnectivityList,1)
    if (TR.Points(TR.ConnectivityList(k,1),3)>0.0 && TR.Points(TR.ConnectivityList(k,2),3)>0.0 && TR.Points(TR.ConnectivityList(k,3),3)>0.0)
    faces = [faces k];
    end
end
TR2.Points = TR.Points;
TR2.ConnectivityList = TR.ConnectivityList(faces,:);
TR2 = triangulation(TR2.ConnectivityList,TR2.Points);
TR2 = SubdivideSphericalMesh(TR2, 1);
TR_final = convhull([TR.Points;TR2.Points]);
tsurf(TR_final,[TR.Points;TR2.Points])
axis equal

F = TR_final;
V = [TR.Points;TR2.Points];

% % OLD CODE (TOO IRREGULAR OUTPUT)
% TR = IcosahedronMesh;
% TR = SubdivideSphericalMesh(TR, 1);
% TR.Points
% TR.ConnectivityList
% P_ico = TR.Points;
% size(P_ico)
% distances = zeros(size(P_ico,1),size(P_handles,1));
% % dsearchn para substituir estes loops
% for i=1:size(P_ico,1)
%     for j = 1: size(P_handles,1)
%         distances(i,j) = norm(P_ico(i,:)-P_handles(j,:));
%     end
% end
% distances = min(distances,[],2);
% probabilities = (1-distances/(max(distances)));
% 
% V = P_ico(rand(size(distances))<=probabilities,:);
% F = convhull(V);
% size(V)
% 
% tsurf(F,V); axis equal; cameratoolbar;
% 
% % [VV,~,FF] = tetgen(V,F,'Flags', '-q1.0/8');
% % figure
% % tsurf(FF,VV); axis equal; cameratoolbar;
% % 
% % 
% % V = VV;
% % F = FF;