classdef bm_cage_edge < handle
    properties
        index
        position
        new_position
    end
    methods
        function obj = bm_cage_edge(index, position)
            obj.index = index;
            obj.position = position;
            obj.new_position = position;
        end
    end
end