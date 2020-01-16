function [new_pos] = pixel2equi(P,imagesize)

new_x = 2*pi*P(:,1)/imagesize(2)-pi;
new_y = -(P(:,2)*pi/imagesize(1)-pi/2);
new_pos = [new_x new_y];
end