function [new_pos] = equi2pixels(P,imagesize)
new_x = (P(:,1)+pi)*imagesize(2)/(2*pi);
new_y = -(P(:,2)-pi/2)*imagesize(1)/pi;
new_pos = [new_x new_y];



end