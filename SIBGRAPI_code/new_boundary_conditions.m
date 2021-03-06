function [b,bc] = new_boundary_conditions(V,F,C,P,E,CE,CB,CCE,bone_list,cage_list,discretization)
  % BOUNDARY_CONDITIONS
  % Compute boundary and boundary conditions for solving for correspondences
  % weights over a set of mesh (V,F) with control points C(p,:) and control
  % bones C(E(:,1),:) --> C(E(:,2),:)
  %
  % [b,bc] = boundary_conditions(V,F,C,P,E,CE)
  %
  % % same as [b,bc] = boundary_conditions(V,F,C,1:size(C,1),[])
  % [b,bc] = boundary_conditions(V,F,C) 
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices (not being used...)
  %  C  list of control vertex positions
  %  P  list of indices into C for point controls, { 1:size(C,1) }
  %  E  list of bones, pairs of indices into C, connecting control vertices, 
  %    { [] }
  %  CE  list of "cage edges", pairs of indices into ***P***, connecting
  %    control ***points***. A "cage edge" just tells point boundary conditions 
  %    to vary linearly along straight lines between the end points, and to be
  %    zero for all other handles. { [] }
  % Outputs
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle ( handle order is point handles then edges handles: P,E)
  % 
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: biharmonic_bounded 
  %

  % control vertices and domain should live in same dimensions
  assert(size(C,2) == size(V,2),'Dims of V and C should match');
  % number o dimensions
  dim = size(C,2);

  % set point handle indices and arrange to be column vector
  if(exist('P','var'))
    if(size(P,1) == 1)
      P = P';
    end
  else
    % if point handles weren't given then treat all control vertices as point
    % handles
    P = (1:size(C,1))';
  end

  % set default edge list to []
  if(~exist('E','var'))
    E = [];
  end

  % set default cage edge list to []
  if(~exist('CE','var'))
    CE = [];
  end

    % set default curved bone list to []
  if(~exist('CB','var'))
    CB = [];
  end
  
  % P should be either empty or a column vector
  assert(isempty(P) || (size(P,2) == 1));
  % E should be empty or be #bones by 2 list of indices
  assert( isempty(E) || (size(E,2) == 2));

  % Modifying here to add C to V
%   V = [V;C];
%   [SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-5);
%   
  % number of point controls
  np = numel(P);
  % number of bone controls
  ne = size(E,1);
  % number of curved bone controls
  ncb = size(CB,1);
  % number of control handles
  m = np + ne + ncb;
  % number of mesh vertices
  n = size(V, 1);

  % number of control vertices
  c = size(C,1);

  if size(C,1) < 100
    %% compute distance from every vertex in the mesh to every control vertex
    % use distances to determine closest mesh vertex to each control vertex
    D = pdist2(V,C);
    % Cv(i) is closest vertex in V to ith control vertex in C
    [minD,Cv] = min(D);
  else
    % Hmmm, this is actually slower for a smaller number of handles 
    Cv = knnsearch(V,C,'K',1);
  end

  % if number of unique closest mesh vertices is less than total number, then
  % we have contradictory boundary conditions
  if(~all(size(unique(Cv)) == size(Cv)))
    warning('Multiple control vertices snapped to the same domain vertex');
  end


  % boundary conditions for all vertices NaN means no boundary conditions
  bc = repmat(NaN,[n m]);

  % compute boundary conditions for point handles
  if( ~isempty(P) )
    bc(Cv(P),:) = eye(np,m);
  end

  sqr_d_tol = 1e-6;

  if isempty(F)
    h = norm(max(V)-min(V))/sqrt(size(V,1));
  else
    h = avgedge(V,F);
  end

  % Compute boundary conditions for bone edge handles
  if(~isempty(E))

    % average edges length to give an idea of "is on edge" tolerance
    
    % loop over bones
    for( ii = 1:size(E,1) )
      % This is projecting to the **sample points**. But if the sample points
      % don't include the end points, this will be totally wrong:
      %[t,sqr_d] = project_to_lines(V,V(Cv(E(ii,1)),:),V(Cv(E(ii,2)),:));
      [t,sqr_d] = project_to_lines(V,C(E(ii,1),:),C(E(ii,2),:));
      on_edge = ((abs(sqr_d) < h*sqr_d_tol) & ((t > -1e-10) & (t < (1+1e-10))));
      % get rid of any NaNs on these rows
      % WARNING: any (erroneous) point handle boundary conditions will get
      % "normalized" with bone boundary conditions
      old = bc(on_edge,:);
      old(isnan(old)) = 0;
      bc(on_edge,:) = old;
      %bc(on_edge,:) = isnan(bc(on_edge,:)).*0 + ~isnan(bc(on_edge,:)).*bc(on_edge,:);
      bc(on_edge,np+ne+ii) = 1;
    end
  end
  
  % Compute boundary conditions for curved bone edge handles
%   if(~isempty(CB))
%     % average edges length to give an idea of "is on edge" tolerance
%         % loop over curved bones
%     for( ii = 1:size(CB,1) )
%       % This is projecting to the **sample points**. But if the sample points
%       % don't include the end points, this will be totally wrong:
%       %[t,sqr_d] = project_to_lines(V,V(Cv(E(ii,1)),:),V(Cv(E(ii,2)),:));
%       %Calcula os coefs do plano que contém o circulo
%       [aplane, bplane, cplane] = bm_planecoeffs(C(CB(ii,1),:),C(CB(ii,2),:));
%       if abs(cplane) <= 1e-5
%         cplane = 1e-4;
%       end
%       %calcula a distancia dos pontos em relacao ao plano
%       [plane_distance] = bm_point2planedist(V,aplane, bplane, cplane, 0);
%       %cria o terceiro ponto que forma o triangulo
%       if C(CB(ii,1),:) == C(CB(ii,2),:)
%          disp('ponto igual') 
%       end
%       [xp, yp, zp] = bm_createtriangle(C(CB(ii,1),:),C(CB(ii,2),:),1);
%       P_trig = [xp yp zp];
%       %projeta os pontos no plano criado
%       V_proj = bm_projpoint2plane(V,aplane,bplane,cplane,0);
%       %calcula as coordenadas baricentricas dos pontos projetados
%       Bar_coords = barycentric_coordinates(V_proj,ones(size(V_proj,1),1)*C(CB(ii,1),:),ones(size(V_proj,1),1)*C(CB(ii,2),:),ones(size(V_proj,1),1)*P_trig);
%       Bar_1 = Bar_coords(:,1);
%       Bar_2 = Bar_coords(:,2);
%       Bar_3 = Bar_coords(:,3);
%       
% %       [t,sqr_d] = project_to_lines(V,C(CB(ii,1),:),C(CB(ii,2),:));
%       t = bm_intersec_retas(V_proj,C(CB(ii,1),:),C(CB(ii,2),:));
%       on_edge = ((t > -1e-5) & (t < 1+1e-5) & (abs(plane_distance) < h) & ...
%     ((Bar_1 > -1e-5) & (Bar_2 > -1e-5) & (Bar_3 > -1e-5) )  & (Bar_1+Bar_2+Bar_3 <= 1+1e-5));
%       %on_edge = ((abs(plane_distance) < h*1e-3) & ((t > -1e-10) & (t < (1+1e-10))));
%        % get rid of any NaNs on these rows
%       % WARNING: any (erroneous) point handle boundary conditions will get
%       % "normalized" with bone boundary conditions
%       old = bc(on_edge,:);
%       old(isnan(old)) = 0;
%       bc(on_edge,:) = old;
%       %bc(on_edge,:) = isnan(bc(on_edge,:)).*0 + ~isnan(bc(on_edge,:)).*bc(on_edge,:);
%       bc(on_edge,np+ii) = 1;
%     end
%   end
  if(~isempty(CB))
    % average edges length to give an idea of "is on edge" tolerance
        % loop over curved bones
    for ii = 1:size(CB,1)
        bone_points = discretization{bone_list(ii)};
        init = bone_points(1:end-1,:);
        fim = bone_points(2:end,:);
        [t,sqr_d] = project_to_lines(V,init,fim);
        [min_dist,ind] = min(sqr_d,[],2);
        tam = size(bone_points,1);
        t = (ind)/tam;
        on_edge = ((abs(min_dist) < h*1e-6) & ((t > -1e-10) & (t < (1+1e-10))));
%         disp([t min_dist on_edge])
%         figure;
%         trisurf(F,V(:,1),V(:,2),V(:,3))
%         hold on
% %         pts = V(on_edge,:)
%         scatter3(V(on_edge,1),V(on_edge,2),V(on_edge,3))
%         plot3(discretization{bone_list(ii)}(:,1),discretization{bone_list(ii)}(:,2),discretization{bone_list(ii)}(:,3))
%         hold off
%       on_edge = dsearchn(V,C(CB(ii,:),:));
%       on_edge = zeros(size(V,1),1);
%       on_edge(indices) = 1;
%       disp(on_edge)

      %on_edge = ((abs(plane_distance) < h*1e-3) & ((t > -1e-10) & (t < (1+1e-10))));
       % get rid of any NaNs on these rows
      % WARNING: any (erroneous) point handle boundary conditions will get
      % "normalized" with bone boundary conditions
      old = bc(on_edge,:);
      old(isnan(old)) = 0;
      bc(on_edge,:) = old;
      %bc(on_edge,:) = isnan(bc(on_edge,:)).*0 + ~isnan(bc(on_edge,:)).*bc(on_edge,:);
      bc(on_edge,np+ii) = 1;
    end
  end

  % compute boundary conditions due to cage edges
  if(~isempty(CE))

    % loop over cage edges
    for( ii = 1:size(CE,1) )
      [t,sqr_d] = project_to_lines(V,C(P(CE(ii,1)),:),C(P(CE(ii,2)),:));
      on_edge = ((abs(sqr_d) < h*sqr_d_tol) & ((t > -1e-10) & (t < (1+1e-10))));
      % get rid of any NaNs on these rows
      % WARNING: clobbers any (erroneous) point handle boundary conditions on
      % points that are on bones)
      old = bc(on_edge,:);
      old(isnan(old)) = 0;
      bc(on_edge,:) = old;
      bc(on_edge,CE(ii,1)) = 1 - t(on_edge);
      bc(on_edge,CE(ii,2)) = t(on_edge);
    end

  end
  
  % compute boundary conditions due to curved cage edges
  if(~isempty(CCE))

    % loop over curved cage edges
    for( ii = 1:size(CCE,1) )
      bone_points = discretization{cage_list(ii)};
        init = bone_points(1:end-1,:);
        fim = bone_points(2:end,:);
        [t,sqr_d] = project_to_lines(V,init,fim);
        [min_dist,ind] = min(sqr_d,[],2);
        tam = size(bone_points,1);
        t = (ind)/tam;
        on_edge = ((abs(min_dist) < h*1e-6) & ((t > -1e-10) & (t < (1+1e-10))));
%         disp(on_edge)
%     figure;
%         trisurf(F,V(:,1),V(:,2),V(:,3))
%         hold on
% %         pts = V(on_edge,:)
%         scatter3(V(on_edge,1),V(on_edge,2),V(on_edge,3))
%         plot3(discretization{cage_list(ii)}(:,1),discretization{cage_list(ii)}(:,2),discretization{cage_list(ii)}(:,3))
%         hold off

%       [t,sqr_d] = project_to_lines(V,C(P(CCE(ii,1)),:),C(P(CCE(ii,2)),:));
%       on_edge = ((abs(sqr_d) < h*sqr_d_tol) & ((t > -1e-10) & (t < (1+1e-10))));
      % get rid of any NaNs on these rows
      % WARNING: clobbers any (erroneous) point handle boundary conditions on
      % points that are on bones)
      old = bc(on_edge,:);
      old(isnan(old)) = 0;
      bc(on_edge,:) = old;
      bc(on_edge,CCE(ii,1)) = 1 - t(on_edge);
      bc(on_edge,CCE(ii,2)) = t(on_edge);
    end

  end

  indices = 1:n;
  % boundary is only those vertices corresponding to rows with at least one non
  % NaN entry
  b = indices(any(~isnan(bc),2));
  bc = bc(b,:);
  % replace NaNs with zeros
  bc(isnan(bc)) = 0;
  % Be sure that any boundary conditions partition unity
  bc(any(bc,2),:)  = bc(any(bc,2),:) ./ repmat(sum(bc(any(bc,2),:),2),1,size(bc,2));
end
