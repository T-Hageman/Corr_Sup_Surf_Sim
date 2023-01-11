function main(caseNum, i)
	meshonly = false;
	maxNumCompThreads(8);

	if (i<=7)
		caseNum = 1;
	else
		caseNum = 2;
		i = i-7;
	end

	switch caseNum
		case 1   % NaCL 1:1000   (i=1:7)
			initO2 = 0.25;
			initNaCL = 10^(1+(i-1)/2);
			initpH = 12;

			Lfrac = 0.5e-2;
			Lx = Lfrac + 1000e-3;
			Hfrac = 2e-2;
			Ly = 0.5e-2;

			OxLim = false;
		case 2   % NaCL 1:1000   (i=8:14 -> i=1:7)
			initO2 = 0.25;
			initNaCL = 10^(1+(i-1)/2);
			initpH = 12;

			Lfrac = 0.5e-2;
			Lx = Lfrac + 1000e-3;
			Hfrac = 2e-2;
			Ly = 0.5e-2;

			OxLim = true;
	end

	sname = "case_" + string(caseNum) + "/i_" + string(i);


	savefolder = "./Results_New/"+sname;
	mkdir(savefolder);
	savefolder=savefolder+"/";
	
	addpath(genpath('./Models'))
	addpath(genpath('./Shapes'))

	delete(gcp('nocreate')); % If no pool, do not create new one.
	p = parpool('threads')
	
	tmr = tic;
	
	%% input properties
	Files=dir(savefolder);
	nmax = 10*(length(Files)-2)-20;
	restartfrom = nmax;

 	if restartfrom>0
 		restart = true;
 		restart_num = restartfrom;
 	else
		restart = false;
		restart_num = 0;
	end

	plotonly = false;
	if restart == false
		% mesh properties
		mesh_in.type = "Fracture";
		mesh_in.dxmin    = 0.1e-3;
		mesh_in.dxmax    = Ly/4;
		mesh_in.Lx    = Lx;
		mesh_in.Ly    = Ly;
		mesh_in.Lfrac = Lfrac;
		mesh_in.Hfrac = Hfrac;
		mesh_in.ipcount1D = 2;
		mesh_in.zeroWeight = true;
		mesh_in.SaveName = savefolder;
		mesh_in.generate = meshonly;
	
		%physics models
		
		physics_in{1}.type = "Electrolyte";
		physics_in{1}.Egroup = "Electrolyte";
		physics_in{1}.D = [9.3; 5.3; 1.3; 2; 1.4; 1; 1]*1e-9;  %H OH Na CL Fe FeOH O2
		physics_in{1}.z = [1; -1; 1; -1; 2; 1; 0];
		physics_in{1}.pH0 = initpH;
		physics_in{1}.NaCl = initNaCL;
		physics_in{1}.O2 = initO2;
		physics_in{1}.Lumped = [true; true]; %water, metal
		physics_in{1}.k = [1e6; 1e-1; 1e-3; 1e-3]; %water, Fe, Fe', FeOH
	
		initH = 1000*10^(-physics_in{1}.pH0);
		initOH = 1000*10^(-14+physics_in{1}.pH0);
		initCl = physics_in{1}.NaCl;
		initNa = initCl-initH+initOH;
	
		F_const = 96485.3329;
		physics_in{2}.type = "ElectrolyteInterface";
		physics_in{2}.Anode = "Anode";
		physics_in{2}.Cathode = "Cathode";
		sf = 1;
		physics_in{2}.k = [	1e-1/F_const*sf,	1e-1/F_const*sf,	0.5,	-0.4;  % Fe <-> Fe2+ (anode)
							1e-4/F_const*sf,	1e-6/F_const*sf,	0.5,	0;     % H+ <->H2    (cathode)
							1e-6/F_const*sf,	1e-6/F_const*sf,	0.5,	0.4;     % O2 <->OH-   (cathode)
							];
		physics_in{2}.ChargeConserve = true;
		physics_in{2}.Em = 0;
		physics_in{2}.Lumped = [1 1 1];
	
		if (OxLim)
			physics_in{3}.type = "OxygenLimiter";
			physics_in{3}.Egroup = "E_Top";
			physics_in{3}.dofs = {"O2"};
			physics_in{3}.conVal = [initO2];
			physics_in{3}.tmax = 48*3600;
		else
			physics_in{3}.type = "Constrainer";
			physics_in{3}.Egroup = "E_Top";
			physics_in{3}.dofs = {"O2"};
			physics_in{3}.conVal = [initO2];
		end

		physics_in{4}.type = "Constrainer";
		physics_in{4}.Egroup = "E_Right";
		physics_in{4}.dofs = {"H";"OH";"Na";"Cl";"Fe";"FeOH"};
		physics_in{4}.conVal = [initH; initOH; initNa; initCl; 0; 0];

		physics_in{5}.type = "Constrainer";
		physics_in{5}.Egroup = "E_Right";
		physics_in{5}.dofs = {"Epot"};
		physics_in{5}.conVal = [0];
	
	
		%% solver inputs
		solver_in.maxIt = 100;
		solver_in.Conv = 1e-4;
		solver_in.tiny = 1e-10;
		solver_in.linesearch = true;
		solver_in.linesearchLims = [0.1 1];
	
		%% initialization
		mesh = Mesh(mesh_in);
		%mesh.plot(false, true, false);
		mesh.check();
		drawnow();

		physics = Physics(mesh, physics_in);
	
	
		dt = 30;
		physics.time = 0;
	
		n_max = 1000*24*360;
		solver = Solver(physics, solver_in);
		tvec = 0;
		I_an_vec = 0;
		I_cat_H_vec = 0;
		I_cat_O_vec = 0;
		Em_vec = 0;
	
		tmax = 60*60*24*14;
	
		startstep = 1;
	else
		filename = savefolder+string(restart_num);
		load(filename, "mesh","physics","solver","dt","tvec","I_an_vec","I_cat_H_vec","I_cat_O_vec","Em_vec","n_max","tmax");
		startstep = restart_num+1;
		solver.linesearch = true;
		solver.linesearchLims = [0.1 1];
		solver.tiny = 1e-10;
		solver.Conv = 1e-4;

		tmax = 60*60*24*100;
	end
	
	if plotonly
		plotres(physics)
	elseif meshonly

	else
		for tstep = startstep:n_max
    		disp("Step: "+string(tstep));
			disp("Time: "+string(physics.time));
			physics.dt = min(3600,dt*1.05^(tstep-1));
			disp("dTime: "+string(physics.dt));
    		
    		solver.Solve();
    		
    		physics.time = physics.time+physics.dt;
    		tvec(end+1) = tvec(end)+physics.dt;
			I_an_vec(end+1) = physics.models{2}.I_anode;
			I_cat_H_vec(end+1)= physics.models{2}.I_Cathode1;
			I_cat_O_vec(end+1)= physics.models{2}.I_Cathode2;
			Em_vec(end+1)   = physics.models{2}.Em;
		
			if mod(tstep, 1) == 0
        		%plotres(physics, tvec, I_an_vec, I_cat_H_vec, I_cat_O_vec, Em_vec);
			end
			if mod(tstep, 10) == 0
        		filename = savefolder+string(tstep);
        		save(filename, "mesh","physics","solver","dt","tvec","I_an_vec","I_cat_H_vec","I_cat_O_vec","Em_vec","n_max","tmax");
			end
	
			if (physics.time>tmax)
				break
			end
		end
		filename = savefolder+"end";
		save(filename, "mesh","physics","solver","dt","tvec","I_an_vec","I_cat_H_vec","I_cat_O_vec","Em_vec","n_max","tmax");
	end
	
	toc(tmr)
end


function plotres(physics, tvec, I_an_vec, I_cat_H_vec, I_cat_O_vec, Em_vec)
    figure(42)
	clf 
    subplot(3,3,1)
        physics.PlotNodal("H",-1, "Electrolyte");
        title("H^+")
		colorbar
	subplot(3,3,2)
        physics.PlotNodal("O2",-1, "Electrolyte");
        title("O_2")
		colorbar
	subplot(3,3,4)
		yyaxis left
		plot(tvec/3600, I_an_vec)
		hold on
		plot(tvec/3600, I_cat_H_vec)
		plot(tvec/3600, I_cat_O_vec)
		legend('Fe','H','O')
		xlabel('t [hours]')
		ylabel('$I_{anode} \;[A]$','Interpreter','latex')
		yyaxis right
		plot(tvec/3600, Em_vec)
		xlabel('t [hours]')
		ylabel('$E_m [V_{SHE}]$','Interpreter','latex')	
	subplot(3,3,3)
        physics.PlotNodal("FeOH",-1, "Electrolyte");
        title("FeOH^{+}")
		colorbar
	subplot(3,3,5)
        physics.PlotNodal("Epot",-1, "Electrolyte");
        title("$\varphi$",'Interpreter','latex')
		colorbar
	subplot(3,3,6)
        physics.PlotNodal("Fe",-1, "Electrolyte");
        title("Fe^{2+}")
		colorbar
	subplot(3,3,7)
        physics.PlotNodal("Na",-1, "Electrolyte");
        title("Na^+")
		colorbar
	subplot(3,3,8)
        physics.PlotNodal("Cl",-1, "Electrolyte");
        title("Cl^-")
		colorbar
	subplot(3,3,9)
        physics.PlotNodal("OH",-1, "Electrolyte");
        title("OH^-")
		colorbar

 	figure(44)
 		clf
 		physics.models{2}.plotReactions(physics);

	figure(43)
		clf
		physics.models{1}.plotFields(physics);

     drawnow();
end

