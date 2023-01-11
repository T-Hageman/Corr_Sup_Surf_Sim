function Solve(obj)

    stop = false;
    it = 0;

	%stepnum = size(obj.convergence_log, 1)+1;
	%cf = figure(9731);
 
    obj.physics.Assemble();
    obj.physics.Constrain();
    recalc_pre=true;
    En_err0 = -1;
    curr_max_it = obj.maxIt;
    while(stop == false)

        fprintf("    Solving it:" + string(it) + "      ");
        tsolve = tic;
        
        recalc_pre = true;
        if (recalc_pre)
            [P,R,C] = equilibrate(obj.physics.K);
            recalc_pre = false;
        end

        if true
            d = -R*P*obj.physics.fint;
            B = R*P*obj.physics.K*C;
			%cond_num(it+1)=condest(B);
			if true
				dy = B\d;
			else
				try
					[L,U] = ilu(B,struct('type','nofill'));
					dy = gmres(B,d,500,1e-4,5000,L,U);
				catch
					dy = B\d;
				end
			end
            dx = C*dy;
        else
            dx = -obj.physics.K\obj.physics.fint;
        end
        tsolve = toc(tsolve);
        fprintf("        (Solver time:"+string(tsolve)+")\n");
		%fprintf("Conditioning numbers: "+string(cond_num(it+1))+"\n");

        if (obj.linesearch && it>0)
            e0 = obj.physics.fint'*dx;
            obj.physics.Update(dx);
            
            obj.physics.Assemble();
            obj.physics.Constrain();
            
            e1 = obj.physics.fint'*dx;
            factor = -e0/(e1-e0);
            factor = max(obj.linesearchLims(1), min(obj.linesearchLims(2), factor));
            obj.physics.Update(-(1-factor)*dx);
            fprintf("    Linesearch: " + string(e0) + " -> " + string(e1) + ":  eta=" + string(factor) +"\n");
        else
            obj.physics.Update(dx);
        end
        
        % convergence
        obj.physics.Assemble();
        obj.physics.Constrain();
        if (En_err0 < 0)
            En_err0 = sum(abs(obj.physics.fint.*dx));
            En_err = En_err0;

            if (En_err0==0)
                En_err0 = 1e-12;
            end
        else
            En_err = sum(abs(obj.physics.fint.*dx));
        end
        En_err_n = En_err/En_err0;

		%obj.convergence_log(stepnum,it+1) = En_err_n;
                    
        fprintf("    Residual:" + string(En_err_n) + "   ("+string(En_err)  +") \n");
        
        it=it+1;
        if (it>curr_max_it || En_err_n<obj.Conv || En_err<obj.tiny)
            obj.physics.Commit("Pathdep");
            irr = obj.physics.Irreversibles();
            if (irr == false)
                stop = true;
            else
                obj.physics.Assemble();
                obj.physics.Constrain();
                recalc_pre = true;
                En_err0 = -1;
                curr_max_it = it + obj.maxIt;
            end
        end
    end
    
    obj.physics.Commit("Timedep");
    %figure(9731);
	%semilogy(obj.convergence_log(stepnum,1:it));
	%hold on
	%drawnow();
end

