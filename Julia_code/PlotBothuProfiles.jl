using JLD2
using DelimitedFiles
using Plots
cd("C:/Users/frank/Documents/GitHub/CAM_DA/Julia_code")
# --- 1. Import Julia Data ---
filename_Julia  = "DAComp1"
filenamet_Julia = filename_Julia * "t.jld2"
filenameU_Julia = filename_Julia * "u.jld2"

@load filenamet_Julia FW_t
@load filenameU_Julia FW_u_RSW

t_Julia = reverse(FW_t)  
u_Julia = FW_u_RSW[end:-1:1, :]

# --- 2. Import Matlab Data ---
filename_Matlab  = "Matlab1"
filenamet_Matlab = filename_Matlab * "_t.csv"
filenameU_Matlab = filename_Matlab * "_u.csv"

# Read the time data and convert it to a 1D vector
FW_t_mat = vec(readdlm(filenamet_Matlab, ',', Float64))

# Read the U data which will automatically parse as an N x 3 Matrix
FW_u_RSW_mat = readdlm(filenameU_Matlab, ',', Float64) 

t_Matlab = FW_t_mat
u_Matlab = FW_u_RSW_mat
u_norms = norm.(eachrow(u_Matlab))

# --- 3. Plotting the Control Profiles ---
num_controls = size(u_Matlab, 2) # Expected to be 3

# Initialize a plot with 3 stacked subplots and linked x-axes
p = plot(layout = (num_controls, 1), size = (800, 700), link = :x, margin=3Plots.mm)
for i in 1:num_controls
    # Plot Julia data (Solid Blue Line)
    xlims!(t_Julia[1],t_Julia[end])
    plot!(p, t_Julia[1:(end-1)], u_Julia[1:(end-1), i], 
          subplot = i, 
          label = (i == 1 ? "Julia" : ""), # Only add legend to the first subplot to prevent clutter
          linewidth = 2, 
          color = :blue,
          title = "Control Profile $i (u_$i)",
          ylabel = "u_$i")
    
          
    # Plot Matlab data (Dashed Red Line)
    plot!(p, t_Matlab, u_Matlab[:, i], 
          subplot = i, 
          label = (i == 1 ? "Matlab" : ""), 
          linewidth = 2, 
          linestyle = :dash, 
          color = :red)
    xlims!(t_Julia[1],t_Julia[end])
end

# Add x-axis label only to the bottom-most subplot
plot!(p, subplot = num_controls, xlabel = "Time [orbits]")

# Display the plot
display(p)