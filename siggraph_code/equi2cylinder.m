function new_position = equi2cylinder(position)
new_position = [cos(position(:,1)) sin(position(:,1)) position(:,2)];
end