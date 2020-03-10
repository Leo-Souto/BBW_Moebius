function [V,F] = adaptive_mesh(P_handles,max_iter)
% ADAPTIVE_MESH Computes a mesh that if finer closer to handles
%
% [V,F] = adaptive_mesh(P_handles,max_iter)
%
% Inputs:
%   P_handles: #P_handles x 3 matrix of handle coordinates
%   max_iter: maximum number of subdivisions of a triangles from the
%   original icosahedral mesh
% Output:
%   V: #V x 3 matrix of 3d vertices
%   F: #F x 3 matrix of indices into V
% 

handle_indices = P_handles(:,2);
P_handles = P_handles(:,1:3);

% 3rd attempt
tic
TR = IcosahedronMesh;
TR = SubdivideSphericalMesh(TR, 2);
disp('initial subdivision')
toc

% tsurf(TR.ConnectivityList,TR.Points)
% axis equal
% hold on
% plot3(P_handles(:,1),P_handles(:,2),P_handles(:,3),'.r','MarkerSize',20);

tic
B = barycenter(TR.Points,TR.ConnectivityList);
disp('barycenter calculation')
toc
% plot3(B(:,1),B(:,2),B(:,3),'.k','MarkerSize',20);
% cameratoolbar
% input('')

tic
% distances = zeros(size(TR.ConnectivityList,1),1);

% for i=1:size(TR.ConnectivityList,1)
%     for j = 1:size(P_handles,1)
% %         distances(i,j) = atan2(norm(cross(B(i,:),P_handles(j,:))),dot(B(i,:),P_handles(j,:)));
%         distances(i,j) = acos(dot(B(i,:),P_handles(j,:)));
% %         distances(i,j) = norm(B(i,:)-P_handles(j,:));
%     end
% end
distances = acos(B*P_handles');
[distances_1,idx] = min(distances,[],2);

distances_2 = zeros(size(idx,1),1);
for k=1:size(idx,1)
    J = find(handle_indices == handle_indices(idx(k)));
    new_distances = distances(k,setdiff(1:size(handle_indices,1),J));
    distances_2(k) = min(new_distances);
end
distances = 0.5*distances_1+0.5*distances_2;
disp('distance calculation')
toc
% figure
% trisurf(TR.ConnectivityList,TR.Points(:,1),TR.Points(:,2),TR.Points(:,3),distances)
% colorbar
% axis equal
% hold on
% plot3(P_handles(:,1),P_handles(:,2),P_handles(:,3),'.r','MarkerSize',20);
% cameratoolbar
% input('')

tic
% subdivisions = round(min(((max_iter/3)*(1./(distances))),max_iter));
subdivisions = round(kumaraswamy(distances,1,max_iter,pi,0));
% subdivisions = round(5*(1-(distances/pi)).^4);
disp('number of extra subdivisions calculation')
toc

tic
for k=1:max_iter
    NewConnectivityList = TR.ConnectivityList(subdivisions==k,:);
    if (~isempty(NewConnectivityList))
        TR_new =  triangulation(NewConnectivityList,TR.Points);
        TR_sub{k} = SubdivideSphericalMesh(TR_new, k);
    else
        TR_sub{k}.Points = [];
    end
end
disp('additional subdivisions')
toc

tic
total_points = 0;
for k=1:max_iter
    total_points = total_points+size(TR_sub{k}.Points,1);
end
disp('total points calculation')
toc

tic
final_points = zeros(total_points,3);
current_number_of_points = 0;
for k=1:max_iter
    final_points(current_number_of_points+1:current_number_of_points+size(TR_sub{k}.Points,1),:) = TR_sub{k}.Points;
    current_number_of_points = current_number_of_points + size(TR_sub{k}.Points,1);
end
disp('selection of final points')
toc

% plot3(final_points(:,1),final_points(:,2),final_points(:,3),'.')
% axis equal

tic
[final_points,~,~] = remove_duplicate_vertices(final_points,1e-7);
disp('remove duplicate')
toc
tic
TR_conv = convhulln(final_points);
disp('convex hull')
toc
% tic
% TR = triangulation(TR_conv,final_points);
% disp('final triangulation')
% toc

F = TR_conv;
V = final_points;

% % 2nd attempt: STILL A BIT IREGULAR
% TR = IcosahedronMesh;
% TR = SubdivideSphericalMesh(TR, 1);
% 
% iter = 1;
% while (iter==1 || (~isempty(faces)) && (iter<3))
%     
%     distances = zeros(size(TR.Points,1),size(P_handles,1));
%     % dsearchn para substituir estes loops
%     for i=1:size(TR.Points,1)
%         for j = 1: size(P_handles,1)
%             distances(i,j) = norm(TR.Points(i,:)-P_handles(j,:));
%         end
%     end
%     distances = min(distances,[],2);
%     
%     faces = [];
%     iter
%     size(TR.ConnectivityList,1)
%     for k=1:size(TR.ConnectivityList,1)
%         if (distances(TR.ConnectivityList(k,1)) + distances(TR.ConnectivityList(k,2)) + distances(TR.ConnectivityList(k,3))<1.5)
%             faces = [faces k];
%         end
%     end
%     faces
%     OldPoints = TR.Points;
%     NewConnectivityList = TR.ConnectivityList(faces,:);
%     if (~isempty(NewConnectivityList))
%         TR2 = triangulation(NewConnectivityList,OldPoints);
%         TR2 = SubdivideSphericalMesh(TR2, 1);
%     end
%     TR_final = convhull([TR.Points;TR2.Points]);
%     TR = triangulation(TR_final,[TR.Points;TR2.Points])
%     iter = iter+1;
% end
% tsurf(TR_final,[TR.Points;TR2.Points])
% axis equal
% cameratoolbar
% 
% F = TR_final;
% V = [TR.Points;TR2.Points];

% % 1st attempt: TOO IRREGULAR OUTPUT
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