classdef bm_handle < handle
    properties
        tag
        position
        new_position
        mesh_position
        new_mesh_position %only used by the curved bones (maybe cages in the future)
        type
        plot_tag
        transformation
        cage_tag %%conjunto de tags que fazem parte de uma mesma cage
        proportion %%proporcao do ponto intermediario ao primeiro
    end
    methods
        function obj = bm_handle(tag, position,mesh_position, type, plot_tag)
            obj.tag = tag;
            obj.position = position;
            obj.new_position = position;
            obj.mesh_position = mesh_position;
            obj.new_mesh_position = mesh_position;
            obj.type = type;
            obj.plot_tag = plot_tag;
        end
    end
end