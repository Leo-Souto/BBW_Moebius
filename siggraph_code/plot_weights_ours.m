weight = pesos.weights(:,:,46);
nx = size(weight,2);
x = linspace(-pi,pi,nx);  
y = linspace(pi/2,-pi/2,nx/2);  
    [X,Y]=meshgrid(x,y);

hold on
surf(X,Y,weight,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5)
colormap(color)
contour3(X,Y,weight,10,'Color','black')
view(0,90);
axis equal
rotate3d off
hold off
axis off
%     axis equal
%     hold on
%     colormap(handles.color)