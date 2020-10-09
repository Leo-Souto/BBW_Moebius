function [V,Vequi,F,index] = bm_createmesh(nx)

ny = floor(nx/2);
dy = cosspace(0,pi,ny);
dy = dy(1:ny/2);
k = size(dy);
dy2 = dy(2:k(2));
dy2 = -dy2;
dy = [fliplr(dy) dy2];
dy = pi/2*dy/max(dy);
NX = floor(nx*cos(dy)) + 2;
points = [];
j = 1;
for a = dy
   for i = linspace(-pi,pi,NX(j))
        points = [points; i a];
   end
   j = j+1;
end
V2 = points;
Vequi = V2;
V = equi2sphere(V2);
[SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-7);
index = knnsearch(SV,V);

F = convhull(SV,'simplify',true);
V=SV;

end