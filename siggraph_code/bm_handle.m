classdef bm_handle < handle
    properties
        tag
        position
        new_position
        mesh_position
        type
        plot_tag
        auxiliar_position
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
            obj.type = type;
            obj.plot_tag = plot_tag;
        end
    end
end