function new_position = cylinder2equi(position)
    new_position = [asin(position(:,2)) position(:,3)];
end