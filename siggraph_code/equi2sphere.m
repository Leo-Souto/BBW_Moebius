function new_position = equi2sphere(position)
    new_position = [cos(position(:,1)).*cos(position(:,2)) sin(position(:,1)).*cos(position(:,2)) sin(position(:,2))];
end