classdef Physics < handle
    %PHYSICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        models
        dofSpace
        
        K
        fint
        StateVec
        StateVec_Old
        time
        dt
        
        condofs
        convals
        convals_corr
        conMat
        unconMat
    end
    
    methods
        PlotNodal(obj, dofName, dispscale, plotloc)
        PlotIP(obj, varName, plotloc)
        
        function obj = Physics(mesh, inputs)
            obj.mesh = mesh;
            obj.dofSpace = DofSpace(obj.mesh);
            
            for i=1:length(inputs)
                f = str2func(inputs{i}.type);
                obj.models{i} = f(mesh, obj, inputs{i});
            end
            
            obj.dt = 0;
            
            dofcount = obj.dofSpace.NDofs;
            obj.StateVec = zeros(dofcount, 1);
            obj.StateVec_Old = obj.StateVec;
            
            obj.K = sparse(dofcount, dofcount);
            obj.fint = zeros(dofcount,1);
        end
        
        function Assemble(obj)
            dofcount = obj.dofSpace.NDofs;

			obj.condofs = [];
			obj.convals = [];

            nonz = round(nnz(obj.K)*1.2);
            obj.K = spalloc(dofcount, dofcount, nonz);
            obj.fint = zeros(dofcount, 1);

            disp("    Assembling:")
            for m=1:length(obj.models)
                obj.models{m}.getKf(obj);
            end
        end
       
        function Commit(obj, commit_type)
            for m=1:length(obj.models)
                obj.models{m}.Commit(obj, commit_type);
            end
            
            if (commit_type == "Timedep")
                obj.StateVec_Old = obj.StateVec;
            end
        end
        
        function anyIrr = Irreversibles(obj)
            anyIrr = false;
            for m=1:length(obj.models)
                anyIrr = anyIrr + obj.models{m}.Irreversibles(obj);
            end
        end
            
        function Constrain(obj)
            obj.convals_corr = obj.convals - obj.StateVec(obj.condofs);
            basemat = speye(size(obj.K));
            obj.unconMat = basemat;
            obj.unconMat(:, obj.condofs) = [];
            obj.conMat = basemat(:, obj.condofs);

            obj.fint = obj.unconMat'*obj.fint - obj.unconMat'*obj.K*obj.conMat*obj.convals_corr;
            obj.K    = obj.unconMat'*obj.K*obj.unconMat;
        end
        
        function Update(obj, dx)
            obj.StateVec = obj.StateVec + obj.unconMat*dx + obj.conMat*obj.convals_corr;
            
        end
        
        function info = Request_Info(obj, var, elems, loc)
            info = false;
            for m=1:length(obj.models)
                [hasInfo, provided] = obj.models{m}.Provide_Info(obj, var, elems, loc);
                if (hasInfo)
                    info = provided;
                end
            end   
        end

    end
end

