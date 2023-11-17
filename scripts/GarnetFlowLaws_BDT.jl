using GarnetFlowLaws
using CairoMakie, MathTeXEngine, GeometryTypes, Enzyme
using Makie.GeometryBasics
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

function main_BDT()

    printfig = true

    # Drucker-Prager
    ϕ  = 35*pi/180
    C  = 5e7
    λf = 0.5

    # P-T-Eii
    Eii = 1e-13
    P   = 1.4e9  # q-coe transition

    P       = LinRange(0.5, 3.5, 400) .* 1e9 
    T_Ji96  = zero(P)
    T_Xu13  = zero(P)
    T_Mei10 = zero(P)
    T_Li06  = zero(P)
    T_Ka95  = zero(P)

    for iP in eachindex(P)
        T_Ji96[iP]  = T_BDT2(Ji96,  Eii, P[iP], λf, ϕ, C) - 273.15
        T_Xu13[iP]  = T_BDT2(Xu13,  Eii, P[iP], λf, ϕ, C) - 273.15
        T_Mei10[iP] = T_BDT2(Mei10, Eii, P[iP], λf, ϕ, C) - 273.15
        T_Li06[iP]  = T_BDT2(Li06,  Eii, P[iP], λf, ϕ, C) - 273.15
        T_Ka95[iP]  = T_BDT2(Ka95,  Eii, P[iP], λf, ϕ, C) - 273.15
    end

    res = 800
    fig = Figure(resolution = (res, res), fontsize=25)
    ax = Axis(fig[1, 1], title = L"$$Garnet brittle-ductile transitions", xlabel = L"$T$ [$^\circ$C]", ylabel = L"$P$ [GPa]", aspect=1.0)
    lines!(ax, T_Ji96,  P, label= "Ji96" , linewidth=4 )
    lines!(ax, T_Xu13,  P, label= "Xu13" , linewidth=4 )
    lines!(ax, T_Mei10, P, label= "Mei10", linewidth=4 )
    lines!(ax, T_Li06,  P, label= "Li06" , linewidth=4 )
    lines!(ax, T_Ka95,  P, label= "Ka95" , linewidth=4 )
    axislegend("References", position = :rt)
    xlims!(ax, 300, 1300)
    display(fig)

    if printfig 
        save( string("figures/Garnet_BDT.png"), fig, px_per_unit = 4) 
    end

    #--------------------------------------------------------------#

end 

main_BDT()
