classdef HydrogenDiffusion < BaseModel
    %Implements lattice hydrogen diffusion through a metal. Input example:
	%	physics_in{2}.type = "HydrogenDiffusion";
	%	physics_in{2}.Egroup = "Metal";
	%	physics_in{2}.DL = 1e-9;
	%	physics_in{2}.NL = 1e6;
    
    properties
        mesh			%pointer to mesh
        myName			%name of this model
        myGroup			%name of element group associated with model
        myGroupIndex	%name of the group index
        dofSpace		%pointer to degree of freedoms
        dofTypeIndices	%vector containing dof numbering
        
        DL			%lattice diffusivity
		NL			%number of lattice sites
        poisson		%poisson ratio
        young		%youngs modulus

		CL_int		%integral of lattice hydrogen
		CL_max		%maximum lattice hydrogen concentration

		R_const = 8.31446261815324;	%gas constant
		T_const = 293.15;			%temperature
		VH_const = 2e-6;			%hydrogen atomic volume
    end
    
    methods
        function obj = HydrogenDiffusion(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "HydrogenDiffusion";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"dx","dy","CL"});
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            %% get parameters
            obj.DL = inputs.DL;
			obj.NL = inputs.NL;
            
			obj.poisson = physics.models{1}.poisson;
            obj.young = physics.models{1}.young;
        end
        
        function getKf(obj, physics)
            fprintf("        HydrogenDiffusion get Matrix:")
            t = tic;
            
            dt = physics.dt;
			Svec = physics.StateVec;
            SvecOld = physics.StateVec_Old;

            dofmatX = [];
            dofmatY = [];
            kmat = [];
			fvec = [];
            dofvec = [];

			CL_sum = 0;		%for integrating hydrogen concentration
			CL_max2 = 0;	%for determining maximum hydrogen concentration
			parfor n_el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
				%get element nodes and shape functions
                Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
				G2 = obj.mesh.getG2(obj.myGroupIndex, n_el);

				% get degrees of freedom and nodal values
                dofsX = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
                dofsY = obj.dofSpace.getDofIndices(obj.dofTypeIndices(2), Elem_Nodes);
                dofsCL= obj.dofSpace.getDofIndices(obj.dofTypeIndices(3), Elem_Nodes);
                dofsXY = [dofsX; dofsY];

                X = Svec(dofsX);
                Y = Svec(dofsY);
                XY = [X;Y];
                CL = Svec(dofsCL);
                CLOld = SvecOld(dofsCL);

				%initialize arrays
                q_el = zeros(length(dofsCL), 1);

                K_cu = zeros(length(dofsCL), length(dofsXY));
                K_cc = zeros(length(dofsCL));
				for ip=1:length(w)
                    %% capacity term
                    q_el = q_el + w(ip) * N(ip,:)'*N(ip,:)*(CL-CLOld)/dt;

                    K_cc = K_cc + w(ip) * N(ip,:)'*N(ip,:)/dt;

                    %% hydraulic stress driven  
					pfx = obj.young/(3*(1-2*obj.poisson));
                    Bstar = obj.getBstar(G2(ip,:,:));
                    dsh = pfx*Bstar*XY;

                   q_el = q_el - w(ip)*obj.DL*obj.VH_const/obj.R_const/obj.T_const * squeeze(G(ip,:,:)) * dsh *max(0,(N(ip,:)*CL));
					K_cc = K_cc - w(ip)*obj.DL*obj.VH_const/obj.R_const/obj.T_const * (squeeze(G(ip,:,:)) * dsh) *N(ip,:);
                   K_cu = K_cu - w(ip)*obj.DL*obj.VH_const/obj.R_const/obj.T_const * squeeze(G(ip,:,:)) * pfx*Bstar*max(0,(N(ip,:)*CL));

                    %% diffusion driven
                    q_el = q_el + w(ip)*obj.DL*(1/max(1e-20,(1-max(0,N(ip,:)*CL)/obj.NL)))*squeeze(G(ip,:,:))*squeeze(G(ip,:,:))'*CL;
                    K_cc = K_cc + w(ip)*obj.DL*(1/max(1e-20,(1-max(0,N(ip,:)*CL)/obj.NL)))*squeeze(G(ip,:,:))*squeeze(G(ip,:,:))' ....
						        + 0*w(ip)*obj.DL*(1/max(1e-20,(1-max(0,N(ip,:)*CL)/obj.NL))^2)*squeeze(G(ip,:,:))*((squeeze(G(ip,:,:))'*CL)*N(ip,:))/obj.NL;

					%% integrations
					CL_sum = CL_sum + w(ip)*N(ip,:)*CL;
					CL_max2 = max(CL_max2, N(ip,:)*CL);
				end

				%% add to sparse allocation vectors
                [dofmatxloc,dofmatyloc] = ndgrid(dofsCL,dofsXY);
                dofmatX = [dofmatX; dofmatxloc(:)];
                dofmatY = [dofmatY; dofmatyloc(:)];
                kmat = [kmat; K_cu(:)];

                [dofmatxloc,dofmatyloc] = ndgrid(dofsCL,dofsCL);
                dofmatX = [dofmatX; dofmatxloc(:)];
                dofmatY = [dofmatY; dofmatyloc(:)];
                kmat = [kmat; K_cc(:)];        

				fvec = [fvec; q_el];
                dofvec = [dofvec; dofsCL];
			end 
			
			% add to stiffness matrix
			physics.fint = physics.fint + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint), 1);
            physics.K = physics.K + sparse(dofmatX, dofmatY, kmat, length(physics.fint),length(physics.fint));
				
			obj.CL_int = CL_sum;
			obj.CL_max = CL_max2;
            
            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
        end
        
		function B = getBstar(~, grad2s) %displacement to hydraulic strain mapping matrix
            cp_count = size(grad2s, 2);
            B = zeros(2, cp_count*2);
            for ii = 1:cp_count %using plane strain e_zz = 0
				%dx
				B(1, ii) = grad2s(1,ii, 1);
				B(2, ii) = grad2s(1,ii, 3);

				%dy
				B(1, ii + cp_count) = grad2s(1,ii, 3);
				B(2, ii + cp_count) = grad2s(1,ii, 2);
            end
		end
        
    end
end

