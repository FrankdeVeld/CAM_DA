using Plots, StaticArrays

function BPlaneTrajectoryPlot(L_Δr_Traj)
    ThPlot          = 1:1:360
    xCirc           = cos.(ThPlot.*pi/180)
    yCirc           = sin.(ThPlot.*pi/180)

    xmin, xmax      = minimum(L_Δr_Traj[:,1])-0.1, maximum(L_Δr_Traj[:,1]) +0.1
    ymin, ymax      = minimum(L_Δr_Traj[:,2])-0.1, maximum(L_Δr_Traj[:,2]) +0.1

    BPlanePlot = Plots.plot(L_Δr_Traj[:,1], L_Δr_Traj[:,2],
        grid = "off",
        legend = false,
        aspect_ratio = :equal,
        linewidth = 2
    )

    Plots.plot!(xCirc, yCirc, # Safe distance
        linewidth = 2
    )
    Plots.xlabel!("b_1 [normalised over σ]")
    Plots.ylabel!("b_2 [normalised over σ]") 
    Plots.title!("Trajectory on the B-plane")
    display(BPlanePlot)
end

function DeltaRMagPlot(L_Δr_Traj, tL)
    normL   = zeros(length(L_Δr_Traj[:,1]))

    for i=1:length(L_Δr_Traj[:,1])
        normL[i] = norm(L_Δr_Traj[i,:])
    end

    NormPlot = Plots.plot(tL, normL,
        grid = "off",
        legend = false,
        linewidth = 2,
    )

    Plots.xlabel!("Time [orbital revolutions]")
    Plots.ylabel!("|Δr| [normalised]") 
    Plots.title!("Time required to reach safety")

    display(NormPlot)
end

function uRSWPlot(tL, u_RSW_L)

    uPlot = Plots.plot(tL[2:(end-1)], u_RSW_L[2:(end-1),1], color=:red, label="R", linewidth=2)
    Plots.plot!(tL[2:(end-1)], u_RSW_L[2:(end-1),2], color=:blue, label="S", linewidth=2)
    Plots.plot!(tL[2:(end-1)], u_RSW_L[2:(end-1),3], grid = "off", color=:green, label="W", linewidth=2)
    
    Plots.xlabel!("t [orbital revolutions]")
    Plots.ylabel!("Thrust magnitude [normalised]") 
    Plots.title!("Control profile throughout time")
    display(uPlot)

    return uPlot
end

function PlotRelDis(t, RelDis, RelVel, TimeDL, PlotBool, JuliaBool)
    row_norms = norm.(eachcol(RelDis))  # Compute row-wise norms
    inner_products = abs.(dot.(eachcol(RelDis), eachcol(RelVel)))
    Startk = 100
    DistanceDL = 7.186789169206707e6
    if PlotBool
        RelDisPlot = Plots.plot(vec(t[end-Startk:end]).*TimeDL, row_norms[end-Startk:end].*DistanceDL, linewidth=3, label="Relative distance")#, yscale=:log10)  

        Plots.xlabel!("Time [s]")
        Plots.ylims!(1500, 2500)
        Plots.ylabel!("Relative distance [m]") 
        if JuliaBool 
            Plots.title!("Relative distance over time, Julia") 
        else 
            Plots.title!("Relative distance over time, Matlab") 
        end

        display(RelDisPlot)

        DotProdPlot = Plots.plot(vec(t[end-Startk:end]).*TimeDL, inner_products[end-Startk:end], linewidth=3, label="Inner product Δr Δv")#, yscale=:log10)  

        Plots.xlabel!("Time [s]")
        Plots.ylabel!("||Inner product|| [-]")  
        if JuliaBool 
            Plots.title!("Inner product Δr Δv over time, Julia")
        else 
            Plots.title!("Inner product Δr Δv over time, Matlab")
        end
        display(DotProdPlot)
    end

    indexmin = argmin(row_norms[end-Startk:end])
    tCA = t[end-Startk+indexmin-1]
    mininner = inner_products[end-Startk+indexmin-1]
    mindis = minimum(row_norms[end-Startk:end])
    return tCA, mindis, mininner
end