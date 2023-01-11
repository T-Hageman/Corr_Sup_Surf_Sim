function [Nodes, Elementgroups, Nodegroups, Area, rect] = Fracture_Generator(obj, props)
    	Lx = props.Lx; 
    	Ly = props.Ly;  
		Lfrac = props.Lfrac;
		Hfrac = props.Hfrac;
    	dxmin = props.dxmin;  
    	dxmax = props.dxmax;   

		sname = props.SaveName;


	if (props.generate)
		%metal
		R1 = [3,4,0,Lx,Lx,0,0,0,Ly,Ly]';
		R2 = [3,4,0,Lfrac,Lfrac,0,0,0,-Hfrac,-Hfrac]';
		gm = [R1,R2];
		sf = '(R1+R2)';
	
		ns = char('R1','R2');
		ns = ns';
		[shp, shpb] = decsg(gm,sf,ns);
		[shp, shpb] = csgdel(shp,shpb,[5]);
	
		figure(754423)
		clf
 		subplot(2,1,1)
 		pdegplot(shp,'EdgeLabels','on','FaceLabels','on','VertexLabels','on')
	
	
		geo = createpde(1);
		geometryFromEdges(geo,shp);
		generateMesh(geo,'Hmax',dxmax,'Hgrad',1.3,'Hedge',{[3,4], dxmin});
	
 		subplot(2,1,2)
 		pdeplot(geo,'NodeLabels','off','ElementLabels','off')
	
		Nodes = geo.Mesh.Nodes';
	
		% interior elements
		Elementgroups{1}.name = "Electrolyte";
		Elementgroups{1}.type = "T6";
		Elementgroups{1}.Elems = geo.Mesh.Elements(:,findElements(geo.Mesh,'region','Face',1))';
	
		Elementgroups{2}.name = "E_Right";
		Elementgroups{2}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[1]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,2));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{2}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
		Elementgroups{3}.name = "E_Top";
		Elementgroups{3}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[2]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,1));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{3}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
		Elementgroups{4}.name = "Cathode";
		Elementgroups{4}.type = "L3B";
		
		cntr = 0;
		N = findNodes(geo.Mesh,'region','Edge',[3]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,2));
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{4}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
		N = findNodes(geo.Mesh,'region','Edge',[5]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,1));
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{4}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
	
	
	
		Elementgroups{5}.name = "Anode";
		Elementgroups{5}.type = "L3B";
		
		N = findNodes(geo.Mesh,'region','Edge',[4]);
		xy = Nodes(N,:,:);
		[~,i] = sort(xy(:,1));
	
		cntr = 0;
		while length(i)>2
			cntr = cntr+1;
			Elementgroups{5}.Elems(cntr,:) = N(i(1:3));
			i(1:2) = [];
		end
    	
    	for g=1:length(Elementgroups)
        	Nodegroups{g}.name = Elementgroups{g}.name;
        	Nodegroups{g}.Nodes = unique(reshape(Elementgroups{g}.Elems,[],1));
		end

		Nodegroups{g+1}.name = "E_CentreTop";
		Nodegroups{g+1}.Nodes= findNodes(geo.Mesh,'region','Vertex',[7]);

		Area = zeros(g,1);
		rect = false;

		save(sname+"Mesh.mat",'Nodes', 'Elementgroups', 'Nodegroups', 'Area', 'rect')
	else
		load(sname+"Mesh.mat",'Nodes', 'Elementgroups', 'Nodegroups', 'Area', 'rect')
	end

end

