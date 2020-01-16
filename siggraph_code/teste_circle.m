nx = 100;
[V,V2,F,indices] = bm_createmesh(nx);
%P1 = equi2sphere(V2(knnsearch(V2,[0 pi/2]),:));
P1 = V(knnsearch(V,equi2sphere([2 0.5])),:);
%P2 = equi2sphere(V2(knnsearch(V2,[0 -pi/2]),:));
P2 = V(knnsearch(V,equi2sphere([0 0])),:);
[aplane, bplane, cplane] = bm_planecoeffs(P1,P2);
if abs(cplane) <= 1e-5
   cplane = 1e-4;
end
[plane_distance] = bm_point2planedist(V,aplane, bplane, cplane, 0)
%       [t,sqr_d] = project_to_lines(V,C(CB(ii,1),:),C(CB(ii,2),:));
V_proj = bm_projpoint2plane(V,aplane,bplane,cplane,0);
[xp, yp, zp] = bm_createtriangle(P1,P2,1);

P_trig = [xp yp zp];

Bar_coords = barycentric_coordinates(V_proj,ones(size(V_proj,1),1)*P1,ones(size(V_proj,1),1)*P2,ones(size(V_proj,1),1)*P_trig);
Bar_1 = Bar_coords(:,1);
Bar_2 = Bar_coords(:,2);
Bar_3 = Bar_coords(:,3);
[t,sqr_d] = project_to_lines(V,P1,P2);

on_edge = ((t > -1e-5) & (t < 1+1e-5) & (abs(plane_distance) < 0.51e-1) & ...
    ((Bar_1 > -1e-5) & (Bar_2 > -1e-5) & (Bar_3 > -1e-5) )  & (Bar_1+Bar_2+Bar_3 <= 1+1e-5));
A = find(on_edge == 1);
sphere
hold on
scatter3(V(A,1),V(A,2),V(A,3),'filled','r')
scatter3([P1(1) P2(1)], [P1(2) P2(2)], [P1(3) P2(3)], 'filled', 'b')
sphere_coord = P1;
sphere_coord2 = P2;
% e_coord = sphere2equi(handles.sphere_points(knnsearch(handles.sphere_points,equi2sphere(equi_coord)),:));
% e_coord2 = sphere2equi(handles.sphere_points(knnsearch(handles.sphere_points,equi2sphere(equi_coord2)),:));
t = linspace(0,1,100);
parametric = [];

parametric(:,1) = (1-t)*sphere_coord(1) + t*sphere_coord2(1);
parametric(:,2) = (1-t)*sphere_coord(2) + t*sphere_coord2(2);
parametric(:,3) = (1-t)*sphere_coord(3) + t*sphere_coord2(3);
norma = sqrt(parametric(:,1).^2 + parametric(:,2).^2 + parametric(:,3).^2);
parametric(:,1) = parametric(:,1)./norma;
parametric(:,2) = parametric(:,2)./norma;
parametric(:,3) = parametric(:,3)./norma;
equi_coord_par = sphere2equi(parametric);
plot3(parametric(:,1),parametric(:,2),parametric(:,3),'LineWidth',2.0)