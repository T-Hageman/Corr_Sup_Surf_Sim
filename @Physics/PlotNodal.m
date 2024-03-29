function PlotNodal(obj, dofName, dispscale, plotloc)
	%Plot nodal data for element group with the name plotloc, for a provided string indicating the dof name, and
	%deforming the plotted results based on the displacemenets with scale
	%dispscale
    dxTypes = obj.dofSpace.getDofType({"dx";"dy";dofName});


    for g=1:length(obj.mesh.Elementgroups)
        if (obj.mesh.Elementgroups{g}.name == plotloc)

            if (obj.mesh.Elementgroups{g}.type == "Q9")
                for el=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                    elnodes = obj.mesh.Elementgroups{g}.Elems(el,:);

                    order = [1 2 3 6 9 8 7 4];

					zdofs = obj.dofSpace.getDofIndices(dxTypes(3), elnodes);
					if (dispscale>0)
                    	xdofs = obj.dofSpace.getDofIndices(dxTypes(1), elnodes);
                    	ydofs = obj.dofSpace.getDofIndices(dxTypes(2), elnodes);

                    	X(el,:) = obj.mesh.Nodes(elnodes(order),1)+dispscale*obj.StateVec(xdofs(order));
                    	Y(el,:) = obj.mesh.Nodes(elnodes(order),2)+dispscale*obj.StateVec(ydofs(order));
					else
                    	X(el,:) = obj.mesh.Nodes(elnodes(order),1);
                    	Y(el,:) = obj.mesh.Nodes(elnodes(order),2);
					end


                    Z(el,:) = obj.StateVec(zdofs(order));
                end
                patch(X'*1000,Y'*1000,Z',Z','EdgeColor','None','FaceColor','interp');
				%fill3(X',Y',Z',Z','FaceColor','interp');

                hold on
                colorbar
			end
			if (obj.mesh.Elementgroups{g}.type == "T6")
                for el=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                    elnodes = obj.mesh.Elementgroups{g}.Elems(el,:);

                    order = [1 4 2 5 3 6];
					order = [1 2 3];

					zdofs = obj.dofSpace.getDofIndices(dxTypes(3), elnodes);
					if (dispscale>0)
                    	xdofs = obj.dofSpace.getDofIndices(dxTypes(1), elnodes);
                    	ydofs = obj.dofSpace.getDofIndices(dxTypes(2), elnodes);

                    	X(el,:) = obj.mesh.Nodes(elnodes(order),1)+dispscale*obj.StateVec(xdofs(order));
                    	Y(el,:) = obj.mesh.Nodes(elnodes(order),2)+dispscale*obj.StateVec(ydofs(order));
					else
                    	X(el,:) = obj.mesh.Nodes(elnodes(order),1);
                    	Y(el,:) = obj.mesh.Nodes(elnodes(order),2);
					end

                    Z(el,:) = obj.StateVec(zdofs(order));
                end
                patch(X'*1000,Y'*1000,Z',Z','EdgeColor','None','FaceColor','interp'); %
				%fill3(X',Y',Z',Z','FaceColor','interp');

                hold on
                colorbar
            end
            if (obj.mesh.Elementgroups{g}.type == "L3" || obj.mesh.Elementgroups{g}.type == "LI6")
                for el=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                    elnodes = obj.mesh.Elementgroups{g}.Elems(el,:);

                    zdofs = obj.dofSpace.getDofIndices(dxTypes(3), elnodes(1:3));

                    order = [1 2 3];
                    X(el,:) = [obj.mesh.Nodes(elnodes(order),1);NaN];
                    Y(el,:) = [obj.mesh.Nodes(elnodes(order),2);NaN];
                    Z(el,:) = [obj.StateVec(zdofs(order));NaN];
                end
                patch(X'*1000,Y'*1000,Z',Z','EdgeColor','interp','FaceColor','None','Marker','o','MarkerFaceColor','flat');
                hold on
                colorbar
			end
			if (obj.mesh.Elementgroups{g}.type == "L3B")
                for el=1:size(obj.mesh.Elementgroups{g}.Elems, 1)
                    elnodes = obj.mesh.Elementgroups{g}.Elems(el,:);

                    zdofs = obj.dofSpace.getDofIndices(dxTypes(3), elnodes(1:3));

                    order = [1 2 3];
                    X(el,:) = [obj.mesh.Nodes(elnodes(order),1);NaN];
                    Y(el,:) = [obj.mesh.Nodes(elnodes(order),2);NaN];
                    Z(el,:) = [obj.StateVec(zdofs(order));NaN];
					Z(el,2) = 0.5*Z(el,2)+0.25*Z(el,1)+0.25*Z(el,3);
                end
                patch(X'*1000,Y'*1000,Z',Z','EdgeColor','interp','FaceColor','None','Marker','o','MarkerFaceColor','flat');
                hold on
                colorbar
            end
        end
    end
end


