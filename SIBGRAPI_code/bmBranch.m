classdef bmBranch < handle
   properties
      isroot = false
      isleaf = false
      filled = false
      xmax
      ymax
      xmin
      ymin
      depth
      leafsID
      leafiterator = 1
      trianglesID
      ntriangles
      branchID
      fatherID
      
   end
   methods 
       %create new leafs
       function [newtree] = bmSpread(obj,tree,V, F)
           obj.isleaf = bmIsLeaf(obj);
           if obj.isleaf %Se for folha vai para o ramo-pai e subdivide ele
               newtree = bmSpread(tree(obj.fatherID),tree,V, F);
           elseif ~obj.isleaf && obj.filled && obj.leafiterator <= 4 %Se não é uma folha e já dividiu passa para sua prox folha
               obj.leafiterator = obj.leafiterator + 1;
               newtree = bmSpread(tree(obj.leafsID(obj.leafiterator - 1)),tree,V, F);
           elseif ~obj.isleaf && obj.leafiterator > 4 && obj.depth > 0 %Passou por todas as folhas, volta para o pai
               newtree = bmSpread(tree(obj.fatherID),tree,V, F);
           elseif obj.leafiterator > 4 && obj.depth == 0 %Voltou tudo e é a raiz, finaliza o loop
               newtree = tree;
           elseif ~obj.isleaf && ~obj.filled && obj.leafiterator <= 4
               leaf1 = bmBranch();
               leaf1.xmax = (obj.xmax-obj.xmin)/2 + obj.xmin;
               leaf1.xmin = obj.xmin;
               leaf1.ymax = obj.ymax;
               leaf1.ymin = (obj.ymax-obj.ymin)/2 + obj.ymin;
               leaf1.depth = obj.depth + 1;
               leaf1.branchID = length(tree) + 1;
               leaf1.fatherID = obj.branchID;
               leaf1.trianglesID = bmSeparateTriangles(obj,leaf1,V, F);
               leaf1.ntriangles = length(leaf1.trianglesID);
               
               
               leaf2 = bmBranch();
               leaf2.xmax = obj.xmax;
               leaf2.xmin = (obj.xmax-obj.xmin)/2 + obj.xmin;
               leaf2.ymax = obj.ymax;
               leaf2.ymin = (obj.ymax-obj.ymin)/2 + obj.ymin;
               leaf2.depth = obj.depth + 1;
               leaf2.branchID = leaf1.branchID + 1;
               leaf2.fatherID = obj.branchID;
               leaf2.trianglesID = bmSeparateTriangles(obj,leaf2,V, F);
               leaf2.ntriangles = length(leaf2.trianglesID);
               
               leaf3 = bmBranch();
               leaf3.xmax = (obj.xmax-obj.xmin)/2 + obj.xmin;
               leaf3.xmin = obj.xmin;
               leaf3.ymax = (obj.ymax-obj.ymin)/2 + obj.ymin;
               leaf3.ymin = obj.ymin;
               leaf3.depth = obj.depth + 1;
               leaf3.branchID = leaf1.branchID + 2;
               leaf3.fatherID = obj.branchID;
               leaf3.trianglesID = bmSeparateTriangles(obj,leaf3,V, F);
               leaf3.ntriangles = length(leaf3.trianglesID);
               
               leaf4 = bmBranch();
               leaf4.xmax = obj.xmax;
               leaf4.xmin = (obj.xmax-obj.xmin)/2 + obj.xmin;
               leaf4.ymax = (obj.ymax-obj.ymin)/2 + obj.ymin;
               leaf4.ymin = obj.ymin;
               leaf4.depth = obj.depth + 1;
               leaf4.branchID = leaf1.branchID + 3;
               leaf4.fatherID = obj.branchID;
               leaf4.trianglesID = bmSeparateTriangles(obj,leaf4,V, F);
               leaf4.ntriangles = length(leaf4.trianglesID);
               
               obj.leafiterator = obj.leafiterator + 1;
               obj.filled = true;
               obj.leafsID = [leaf1.branchID leaf2.branchID leaf3.branchID leaf4.branchID];
               tree = [tree leaf1 leaf2 leaf3 leaf4];
               newtree = bmSpread(tree(leaf1.branchID),tree,V, F);
              
           end
       end
       
       %Leaf Test
       function value = bmIsLeaf(obj)
           if obj.isleaf
              value = true;
           end
           if ~obj.isleaf
               if obj.depth > 5 || obj.ntriangles < 5
                   value = true;
               else
                   value = false;   
               end
           end
       end
   end
end