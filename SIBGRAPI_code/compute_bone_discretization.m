function [Position, t_value] = compute_bone_discretization(handle,Transf)
% This function returna the position and proportion of a
% discretization of a bone handle

a = Transf(1);
b = Transf(2);
c = Transf(3);
d = Transf(4);

spos = equi2sphere(handle.position);
t_value = linspace(0,1,50)';
omega = acos(dot(spos(1,:),spos(3,:)));
spos2x = sin((1-t_value)*omega)/sin(omega)*spos(1,1) + sin((t_value)*omega)/sin(omega)*spos(3,1);
spos2y = sin((1-t_value)*omega)/sin(omega)*spos(1,2) + sin((t_value)*omega)/sin(omega)*spos(3,2);
spos2z = sin((1-t_value)*omega)/sin(omega)*spos(1,3) + sin((t_value)*omega)/sin(omega)*spos(3,3);
spos2 = [spos2x spos2y spos2z];

% spos2 = slerp(spos(1,:),spos(3,:),t_value);
X_C =(2*spos2(:,2))./(spos2(:,1) + 1);
Y_C =(2*spos2(:,3))./(spos2(:,1) + 1);
Complex = X_C + Y_C*1i;
Complex2 = (a*Complex + b)./(c*Complex + d);
X_C_n = real(Complex2);
Y_C_n = imag(Complex2);
    
r_polar = sqrt(X_C_n.^2 + Y_C_n.^2);
theta_polar = atan2(X_C_n,Y_C_n);
theta_polar = theta_polar + pi/2;
X_C_n = r_polar.*cos(theta_polar);
Y_C_n = r_polar.*sin(theta_polar);
    
x_sphere2 = -(X_C_n.^2 + Y_C_n.^2 - 4)./(X_C_n.^2 + Y_C_n.^2 + 4);
y_sphere2 = -(4*X_C_n)./(X_C_n.^2 + Y_C_n.^2 + 4);
z_sphere2 = (4*Y_C_n)./(X_C_n.^2 + Y_C_n.^2 + 4);
sphere_final = [x_sphere2  y_sphere2 z_sphere2];
Position = sphere2equi(sphere_final);

end