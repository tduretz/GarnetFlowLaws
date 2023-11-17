using GarnetFlowLaws
using CairoMakie, MathTeXEngine, GeometryTypes, Enzyme
using Makie.GeometryBasics
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

function main_Zulauf2023()

    printfig = true

    # Drucker-Prager
    ϕ  = 35*pi/180
    C  = 5e7
    λf = 0.5
    
    #--------------------------------------------------------------#
    # Reproduce figure from Zulauf et al 202+
    e      = 1e-13
    T_ax   = LinRange(100, 500, 400) .+ 273.15 
    P_ax   = LinRange(0.0, 1.9, 400) .* 1e9 
    fH2O   = zeros(length(T_ax), length(P_ax))
    T_Wu13 = zeros(length(P_ax))
    
    # Compute water fugracity for each T and P
    for iP in eachindex(P_ax)
        for iT in eachindex(T_ax)
            fH2O[iT,iP]  = PSfugacity(P_ax[iP], T_ax[iT])
        end
    end 
    fH2O[isnan.(fH2O)] .= 1.
    fH2O[isinf.(fH2O)] .= 1.

    @show minimum(fH2O)
    @show maximum(fH2O)

    σ_Wu13 = σ(Xu13,  e, T_ax, P_ax', fH2O)

    σ_brit = 2*(C.*cos(ϕ) .+ P_ax'.*sin(ϕ)*λf .+ 0 .* T_ax)

    σ_Wu13[σ_Wu13.>σ_brit] .= σ_brit[σ_Wu13.>σ_brit]
    for iP=1:length(P_ax)
        T_Wu13[iP] = T_BDT2(Xu13,  e, P_ax[iP], λf, ϕ, C) .- 273.15
    end

    # box 
    Pmin = 1.0
    Pmax = 1.2
    Tmin = 320.
    Tmax = 400.
    Tr   = LinRange(Tmin, Tmax, 100)
    Pr   = LinRange(Pmin, Pmax, 100)

    res = 800
    fig = Figure(resolution = (res, res), fontsize=25)
    ax = Axis(fig[1, 1], title = L"$$Garnet deformation map, $\dot{\varepsilon} = 10^{-13}$ s$^{-1}$", xlabel = L"$T$ [°C]", ylabel = L"$P$ [GPa]", aspect=1.0)
    hm=heatmap!(ax, T_ax.-273.15, P_ax./1e9, log10.(σ_Wu13) , colormap=cgrad(:roma, rev = true))
    contour!(ax, T_ax.-273.15, P_ax./1e9, σ_Wu13, levels=[700e6, 500e6, 300e6, 100e6], color=:white)
    lines!(ax, T_Wu13, P_ax./1e9, color=:white, linewidth=10 )
    lines!(ax, Tr, Pmin.*ones(size(Tr)), color=:blue, linewidth=4, linestyle=:dash )
    lines!(ax, Tr, Pmax.*ones(size(Tr)), color=:blue, linewidth=4, linestyle=:dash )
    lines!(ax, Tmin.*ones(size(Pr)), Pr, color=:blue, linewidth=4, linestyle=:dash )
    lines!(ax, Tmax.*ones(size(Pr)), Pr, color=:blue, linewidth=4, linestyle=:dash )
    Colorbar(fig[1,2], hm, width = 20, label = L"$\log_{10} (\sigma_1 - \sigma_3)$ [Pa]",
    labelsize = 25, ticklabelsize = 20, bbox=ax.scene.px_area,
    alignmode = Outside(10), halign = :left, ticklabelcolor = :black, labelcolor = :black,
    tickcolor = :white)

    text!(ax, 450, 0.5, text=L"$100$ MPa", rotation=π/2-0.13, color=:white)
    text!(ax, 396, 0.70, text=L"$300$ MPa", rotation=π/2-0.13, color=:white)
    text!(ax, 385, 1.21, text=L"$500$ MPa", rotation=π/2-0.13, color=:white)
    text!(ax, 371, 1.23, text=L"$700$ MPa", rotation=π/2-0.13, color=:white)

    # Frictional box
    xbox = [215., 317, 317., 215]
    ybox = [0.15, 0.15, 0.35, 0.35] .+ 1.1
    p = [Polygon( Point2f[ (xbox[j], ybox[j]) for j=1:4] ) for i in 1:1]
    poly!(ax, p, color = :white)
    text!(ax, 228, 1.36, text=L"$$Frictional", color=:black, fontsize=32)
    text!(ax, 235, 1.27, text=L"$\lambda_\mathrm{fluid} = 0.5$", color=:black, fontsize=28)
    
    # Viscous box
    xbox = [400., 490, 490., 400]
    ybox = [0.15, 0.15, 0.35, 0.35] .+ 0.05
    p = [Polygon( Point2f[ (xbox[j], ybox[j]) for j=1:4] ) for i in 1:1]
    poly!(ax, p, color = :white)
    text!(ax, 416, 0.3, text=L"$$Viscous", color=:black, fontsize=32, rotation=0(π/2-0.13))
    text!(ax, 407, 0.21, text=L"$f_\mathrm{H2O} = f(T,P)$", color=:black, fontsize=28, rotation=0(π/2-0.13))

    xlims!(ax, 200., 500. )
    ylims!(ax, 0., 1.5 )
    if printfig 
        save( string("figures/GarnetDeformationMap_Xu13_lamf$(λf).png"), fig, px_per_unit = 4) 
    end
    display(fig)
    #--------------------------------------------------------------#
end 

main_Zulauf2023()
