classdef DofSpace < handle
    %DOFSPACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        DofTypes
        DofNumbering   %node, dofType
        NDofs
		Scales
		ScalesF
    end
    
    methods
        function obj = DofSpace(mesh)
            %DOFSPACE Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.mesh = mesh;
            obj.DofTypes = {};
            obj.NDofs = 0;
			obj.Scales = [];
			obj.ScalesF = [];
            obj.DofNumbering = sparse(length(obj.mesh.Nodes), length(obj.DofTypes));
        end
        
        function dofIndex = addDofType(obj, dofnames)
            dofIndex = zeros(length(dofnames),1);
            dofTypeCount = length(obj.DofTypes);
            
            for i=1:length(dofnames)
               alreadyexist = false;
               for j=1:dofTypeCount
                  if (dofnames{i}==obj.DofTypes{j})
                      alreadyexist = true;
                      dofIndex(i) = j;
                  end
               end
               if (alreadyexist == false)
                   dofTypeCount = dofTypeCount+1;
                   obj.DofTypes{dofTypeCount} = dofnames{i};
                   dofIndex(i) = dofTypeCount;
                   obj.DofNumbering(:,length(obj.DofTypes)) = 0;
				   obj.Scales = [obj.Scales; 1];
				   obj.ScalesF = [obj.ScalesF; 1];
               end
            end
            
            
        end
        

        function addDofs(obj, dofIndices, nodeIndex)
            for i=1:length(dofIndices)
                for k=1:length(nodeIndex)
                    if (nodeIndex(k)<=size(obj.DofNumbering, 1))
                        curNum = obj.DofNumbering(nodeIndex(k), dofIndices(i));
                    else
                        obj.DofNumbering(nodeIndex(k), :) = 0;
                        curNum = 0;
                    end
                    if (curNum == 0) %dof does not yet exist if 0, add dof
                        obj.NDofs = obj.NDofs+1;
                        obj.DofNumbering(nodeIndex(k), dofIndices(i)) = obj.NDofs;
                    end
                end
            end
        end
        
        function DofTypeIndex = getDofType(obj, dofnames)
            DofTypeIndex = zeros(length(dofnames),1);
            for i=1:length(dofnames)
               for j=1:length(obj.DofTypes)
                  if (dofnames{i}==obj.DofTypes{j})
                      DofTypeIndex(i) = j;
                  end
               end
            end
        end
        
        function DofIndices = getDofIndices(obj, dofType, NodeIndices)
           
           DofIndices = obj.DofNumbering(NodeIndices, dofType);
           if (length(find(DofIndices==0)) ~= 0)
              disp("error here (dofspace)") 
           end
		end

		function scale = getScale(obj, dofType)
			scale = obj.Scales(dofType);
		end

		function setScale(obj, dofnames, value)
			typenums = obj.getDofType(dofnames);
			obj.Scales(typenums) = value;
		end

		function scale = getScaleF(obj, dofType)
			scale = obj.ScalesF(dofType);
		end

		function setScaleF(obj, dofnames, value)
			typenums = obj.getDofType(dofnames);
			obj.ScalesF(typenums) = value;
		end

    end
end

