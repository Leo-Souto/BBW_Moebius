function new_position = sphere2equi(position)
new_position = [atan2(position(:,2), position(:,1)) asin(position(:,3))];
end