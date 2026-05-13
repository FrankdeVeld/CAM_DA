function input = write_input(scenario, params)

    input         = struct();
    input.N       = params.nx_orb*params.n_orb+1; 
    input.Lsc     = scenario.Lsc; 
    input.et      = 659871.07119168108; 
    input.t_back  = params.n_orb*2*pi; 
    input.uMax    = scenario.ctrlMax;
    input.xp_tCA  = scenario.x_p';
    input.xs_tCA  = scenario.x_s';
    input.P       = scenario.cov;
    input.HBR     = scenario.HBR;
    input.metric_case = 2;
    input.breakOnThreshold = 0;
    if input.metric_case == 1
        input.lim = (scenario.md_lim/Lsc)^2;
    else
        input.lim = scenario.smdLim;
    end
    input.tCAHandling = 2;
    
    fid = fopen('./input.json','w'); 
    fwrite(fid,jsonencode(input),'char'); 
    fclose(fid);

end