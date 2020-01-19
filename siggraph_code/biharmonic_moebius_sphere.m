function [new_image] = biharmonic_moebius_sphere(Image,V,F,handles, T_coefs, size_final)
%HELP Written here
%% The mesh construction takes place before the image editing
% We assume an already subdivided icosahedral mesh as input
% Lets assume only point handles at this time, we must pass all the handles
% data to the function
num_of_handles = length(T_coefs(:,1));
positions = [];
is_conform = [];

for i = 1:num_of_handles
    switch handles{i}.type
        case 'Point'
            positions = cat(1,positions,handles{i}.new_position(1,:));
        case 'Curved'
            positions = cat(1,positions,handles{i}.new_mesh_position);
    end
end
point_indices = [];
curved_bone_indices = [];
for i = 1:num_of_handles
    switch handles{i}.type
        case 'Point'
            indice = dsearchn(positions,handles{i}.new_position(1,:));
            point_indices = [point_indices,indice];
%             new_handle_order{end+1} = handles{i};
        case 'Curved'
            indice = dsearchn(positions,handles{i}.new_mesh_position);
            curved_bone_indices = [curved_bone_indices; indice'];
    end
end
V = [V; equi2sphere(positions)];
F = convhull(V);
%computing spherical boundary conditions and biharmonic weights
[b,bc] = new_boundary_conditions(V,F,equi2sphere(positions),point_indices,...
    [],[],curved_bone_indices,[]);
if ~isempty(is_conform)
    W = biharmonic_bounded(V,F,b,bc,'POU',false,...
        'ShapePreserving',is_conform);
else
    W = biharmonic_bounded(V,F,b,bc,'POU',false);
end
soma = sum(W,2);
for i=1:num_of_handles
   W(:,i) = W(:,i)./soma; 
end
% figure;
% trisurf(F,V(:,1),V(:,2),V(:,3),W(:,4))
%% irregular mesh to grid weight interpolation

V_equi = sphere2equi(V);
size(V_equi)
[SV,SVI,~] = remove_duplicate_vertices(V_equi,1e-7);
size(SV)
size(W)
index = dsearchn(SV,V_equi);
W = W(SVI,:);
size(W)
nx = 150;
x = linspace(-pi+1e-6,pi-1e-6,nx);
y = linspace(pi/2-1e-6,-pi/2+1e-6,nx/2);
[X,Y]=meshgrid(x,y);
%%%%%%%%tem que trocar griddata por scattered interpolant pra ganhar
%%%%%%%%performance
Interpolant = scatteredInterpolant(SV,W(:,1));
non_expanded = zeros(nx/2,nx,num_of_handles); %from irregular to regular mesh
for i  = 1:num_of_handles
    Interpolant.Values = W(:,i);
    non_expanded(:,:,i) = Interpolant(X,Y);
end
% h = figure;
% figure(h)
% trisurf(F,V(:,1),V(:,2),V(:,3),W(:,1))
%from small regular, to final size weights
nx = size_final(2);
ny = size_final(1);
x = linspace(-pi+1e-6,pi-1e-6,nx);
y = linspace(pi/2-1e-6,-pi/2+1e-6,ny);
[U,V]=meshgrid(x,y);
W_expanded = zeros(ny,nx,num_of_handles); %expand to image size
for i  = 1:num_of_handles
    W_expanded(:,:,i) = qinterp2(X,Y,non_expanded(:,:,i),U,V);
end
%%%%%finished weights computation%%%%%

%% Image editing
%desconcatena a imagem
R = double(Image(:,:,1));
G = double(Image(:,:,2));
B = double(Image(:,:,3));

%pega o tamanho da imagem
[nr ,nc]=size(R);

%Define the domain
aa = linspace(-pi,pi,nc);
cc = linspace(-pi/2,pi/2,nr);
% [AA,CC] = meshgrid(aa,cc);
[AA,CC]=ndgrid(aa,cc);

%Define the counter-domain
x = linspace(-pi,pi,size_final(1,2));  
y = linspace(pi/2,-pi/2,size_final(1,1));
% x = V(:,1);
% y = V(:,2);
[X,Y]=ndgrid(x,y);

X_total = zeros(size_final(1,2),size_final(1,1));
Y_total = zeros(size_final(1,2),size_final(1,1));
Z_total = zeros(size_final(1,2),size_final(1,1));


%from final equirectangular projection to the sphere
x_sphere = cos(X).*cos(Y);
y_sphere = sin(X).*cos(Y);
z_sphere = sin(Y);
%From sphere to stereographic plane
X_C =(2*y_sphere)./(x_sphere + 1);
Y_C =(2*z_sphere)./(x_sphere + 1);
%From stereographic plane to complex numbers
Complex = X_C + Y_C*1i;


for i = 1:num_of_handles
    Current_Weight = W_expanded(:,:,i)';
    a = T_coefs(i,1);
    b = T_coefs(i,2);
    c = T_coefs(i,3);
    d = T_coefs(i,4);
    if abs(c) < 1e-10
        c = 1e-10 + 1e-10*1i;
    end
    %inverse moebius transformation
    Complex2 = Complex - a/c;
    Complex2 = Complex2*(-(c^2/(a*d - b*c)));
    Complex2 = 1./Complex2;
    Complex2 = Complex2 - d/c;
    %Returning to stereographic plane
    X_C_n = real(Complex2);
    Y_C_n = imag(Complex2);
    %Fixing the mirrowing
    r_polar = sqrt(X_C_n.^2 + Y_C_n.^2);
    theta_polar = atan2(X_C_n,Y_C_n);
    theta_polar = theta_polar + pi/2;
    X_C_n = r_polar.*cos(theta_polar);
    Y_C_n = r_polar.*sin(theta_polar);
 
    %Return to the sphere
    x_sphere2 = -(X_C_n.^2 + Y_C_n.^2 - 4)./(X_C_n.^2 + Y_C_n.^2 + 4); 
    y_sphere2 = -(4*X_C_n)./(X_C_n.^2 + Y_C_n.^2 + 4); 
    z_sphere2 = (4*Y_C_n)./(X_C_n.^2 + Y_C_n.^2 + 4);
    
    Delta_X = Current_Weight.*x_sphere2;
    Delta_Y = Current_Weight.*y_sphere2;
    Delta_Z = Current_Weight.*z_sphere2;
%     if i == 2 
%         disp([Delta_X])
%     end
    X_total = X_total + Delta_X;
    Y_total = Y_total + Delta_Y;
    Z_total = Z_total + Delta_Z;
   
end


U_final = real(atan2(Y_total, X_total));
V_final = real(-asin(Z_total));

R = R';
G = G';
B = B';

R_interp = griddedInterpolant(AA,CC,R);
G_interp = griddedInterpolant(AA,CC,G);
B_interp = griddedInterpolant(AA,CC,B);

r = R_interp(U_final,V_final)';
g = G_interp(U_final,V_final)';
b = B_interp(U_final,V_final)';

%final image (concatenating matrix layers)
new_image = (uint8(cat(3,r,g,b)));

end