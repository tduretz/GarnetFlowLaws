using CairoMakie, MathTeXEngine, GeometryTypes
Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

# using SymPy
# @vars R n T P E V B C ϕ ε̇ σ f r
# F = ε̇ - A*σ^n*f^r*exp(- (E + P*V)/R/T)
# σv = solve(F, σ)
# display(σv)
printfig = false

# Gas constant
R   = 8.316

# Shear modulus
G  = 8e10
Tm = 1853/0.93

# Drucker-Prager
phi = 35*pi/180
Coh  = 1e7
fO2  = 1e-17

d = 2.5e-3

T_BDT(p, Eii, P, fH2O) =  (p.Q .+ P .* p.V) ./ (R .* log.(p.A .* fH2O .^ p.r .* (Coh .* cos(phi) .+ P .* sin(phi)) .^ p.n ./ Eii))

ε̇(p, σ, T, P, d, fO2)  = 1.0./(T .* d.^p.m) .* p.A*σ.^p.n*fO2^p.r*exp(- (p.Q + P*p.V)/R/T)

σdiff(p, ε̇, T, P, d, fO2) = (T .* d.^p.m .* ε̇ .* fO2 .^ (-p.r) .* exp.((p.Q .+ P .* p.V) ./ (R .* T)) ./ p.A) .^ (1 ./ p.n)

σdis(p, ε̇, T, P, fH2O) = (ε̇ .* fH2O .^ (-p.r) .* exp.((p.Q .+ P .* p.V) ./ (R .* T)) ./ p.A) .^ (1 ./ p.n)


function main()

# P-T-Eii
T   = 688+273
Eii = 1e-13
P   = 2.8e9  # q-coe transition

# Creep - Ji & Matignole (1996)
Ji96 = (
    n    = 2.22,
    Q    = 485e3,
    A    = 5.84e6 *10^(-6*2.22), # from MPa^(-n)/s to Pa^(-n)/s ---------> 2.7952e-07 --- Gut!
    r    = 0.0,
    V    = 0.0
)

# Creep - Xu et al. (2013)
Xu13 = (
    n = 3.5,
    r = 1.0,
    Q = 215e3,
    A = 10^(-5.3) *10^(-6*(3.5 + 1.0)), # from MPa^(-n)/s to Pa^(-n)/s
    V = 28e-6
)

# Diffusion creep - Wang et al. (2010)
Wang00 = (
    n = 1.0,
    m = 2.5,
    Q = 347e3,
    A = 5.32e-6*10^(-6*1.1),
    r = 0.0,
    V = 0.0
)


#--------------------------------------------------------------#
fH2O = 1.2e9
# Reproduce figure from Zulauf et al 202+
e      = 1e-13
d_ax   = 10 .^LinRange(-4, -1, 400)  
T_ax   = LinRange(100, 500, 400) .+ 273.15 
P_ax   = LinRange(0.0, 1.5, 400) .* 1e9 
σ_Xu13 = σdis(Xu13,  e, T_ax, P_ax', fH2O)
# σ_Wa00 = σdiff(Wang00, Eii, T_ax, P_ax', d, fO2)
σ_brit = Coh.*cos(phi) .+ P_ax'.*sin(phi) .+ 0 .* T_ax

# σ_Wa00 = 0*P_ax' .+ (e .* (5.32e-6*10^(-6*1.1)) ./T_ax * (2.5e-3)^(-2.5) .*exp.(-347000 ./8.314 ./T_ax ) ).^(-1/1.1) 

# σ_Wa00      = 1e6.*ones(length(T_ax), length(P_ax))
# for i =1:100
#     E_Wa00   = 0*P_ax' .+ (σ_Wa00.^1.1 .* (5.32e-6) ./T_ax * (2.5e-3)^(-2.5) .*exp.(-347000 ./8.314 ./T_ax ) )
#     f        = E_Wa00 .- 1e-13
#     σ_Wa00 .-= f.*1e17
# end

# σ = 1e8
# E_Wa00   = 0*P_ax' .+ (σ.^1.1 .* (5.32e-6) ./T_ax * (2.5e-3)^(-2.5) .*exp.(-347000 ./8.314 ./T_ax ) )
# @show E_Wa00[end]
# sss = (E_Wa00[end] .* (5.32e-6)^(-1) .* T_ax * (2.5e-3)^(2.5) .*exp.(347000 ./8.314 ./T_ax ) ).^(1/1.1)
# @show sss[end]
# @show σ-sss[end]

d = 100e-6
σ_Wa00 = 0*P_ax' .+ (e .* (5.32e-6)^(-1) .* T_ax .* (d)^(2.5) .*exp.(347000 ./8.314 ./T_ax ) ).^(1/1.1) 

@show size(σ_Wa00)

dfm = zero(σ_brit)
dfm[σ_Xu13.<σ_Wa00 .&& σ_Xu13.<σ_brit] .= 1
dfm[σ_Xu13.>σ_Wa00 .&& σ_Wa00.<σ_brit] .= 2
@show sum(dfm.==0)
@show sum(dfm.==1)
@show sum(dfm.==2)


σ_Xu13[σ_Xu13.>σ_brit] .= σ_brit[σ_Xu13.>σ_brit]
σ_Xu13[σ_Xu13.>σ_Wa00] .= σ_Wa00[σ_Xu13.>σ_Wa00]
T_Xu13 = T_BDT(Xu13,  e, P_ax, fH2O) .- 273.15

 # box 
 Pmin = 1.0
 Pmax = 1.2
 Tmin = 320.
 Tmax = 400.
 Tr   = LinRange(Tmin, Tmax, 100)
 Pr   = LinRange(Pmin, Pmax, 100)

res = 800
fig = Figure(resolution = (res, res), fontsize=25)
ax = Axis(fig[1, 1], title = L"$$Garnet deformation map (Wang et al., 2000), $d = 100$ μm", xlabel = L"$T$ [°C]", ylabel = L"$P$ [GPa]", aspect=1.0)
heatmap!(ax, T_ax.-273.15, P_ax./1e9, log10.(σ_Wa00) , colormap=(:lajolla))
# heatmap!(ax, T_ax.-273.15, P_ax./1e9,dfm , colormap=(:lajolla))

# heatmap!(ax, T_ax.-273.15, P_ax./1e9, log10.(E_Wa00) , colormap=(:lajolla))
# contour!(ax, T_ax.-273.15, P_ax./1e9, E_Wa00, levels=[1e-14, 1e-13], color=:white)


contour!(ax, T_ax.-273.15, P_ax./1e9, σ_Wa00, levels=[800e6, 400e6, 200e6, 100e6], color=:black)

@show minimum(σ_Wa00)
@show maximum(σ_Wa00)

lines!(ax, Tr, Pmin.*ones(size(Tr)), color=:blue, linewidth=3, linestyle=:dash )
lines!(ax, Tr, Pmax.*ones(size(Tr)), color=:blue, linewidth=3, linestyle=:dash )
lines!(ax, Tmin.*ones(size(Pr)), Pr, color=:blue, linewidth=3, linestyle=:dash )
lines!(ax, Tmax.*ones(size(Pr)), Pr, color=:blue, linewidth=3, linestyle=:dash )

# lines!(ax, T_Xu13, P_ax./1e9, color=:white, linewidth=10 )
# poly!(ax, [Rectangle{Float64}(125, 125, 50, 50)], color = :gray, strokewidth = 2, strokecolor = :red)
# text!(ax, 300, 0.7, text=L"$$Brittle-Ductile transition", rotation=π/2-0.15, color=:white)
# text!(ax, 400, 0.7, text=L"$100$ MPa", rotation=π/2-0.33, color=:white)
# text!(ax, 380, 0.9, text=L"$200$ MPa", rotation=π/2-0.33, color=:white)
# text!(ax, 375, 1.3, text=L"$400$ MPa", rotation=π/2-0.33, color=:white)
# text!(ax, 365, 1.6, text=L"$800$ MPa", rotation=π/2-0.30, color=:white)

# xlims!(ax, 100., 500. )
# ylims!(ax, 0., 1.5 )
if printfig save( string(@__DIR__,  "/GarnetDeformationMap_Xu13.png"), fig, px_per_unit = 4) end
display(fig)




# Creep pyrope - Li et al. 2006

# Flow stress
# f = 1.0#1.0/6.0*2.0^(1.0/n) * 3.0^((n-1.0)/2.0/n)
# B = f*A^(-1/n) * exp(Q/n/R/T)
# C = (2*B)^-n
# S = B*Eii^(1/n)

# Find BDT
# @show Tbdt = Q./(R.*n.*log(A.^(1.0./n).*Eii.^(-1.0./n).*(Coh.*cos(phi) + P.*sin(phi))./f)) - 273.15

# fH2O = 12e9

# @show T_BDT(Ji96,  Eii, P, fH2O) - 273.15
# @show T_BDT(Xu13,  Eii, P, fH2O) - 273.15
# @show T_BDT(Mei10, Eii, P, fH2O) - 273.15
# @show T_BDT(Li06,  Eii, P, fH2O) - 273.15
# @show T_BDT(Ka95,  Eii, P, fH2O) - 273.15

# #--------------------------------------------------------------#
# # Reproduce figure 8 from Xu et al 2013
# logσ    = 1:0.1:3
# σ_ax    = 10.0.^logσ
# ε̇_Xu13  = ε̇(Xu13,  σ_ax.*1e6, 1473., 2e9, fH2O)
# ε̇_Mei10 = ε̇(Mei10, σ_ax.*1e6, 1473., 2e9, fH2O)
# ε̇_Li06  = ε̇(Li06,  σ_ax.*1e6, 1473., 2e9, fH2O)
# ε̇_Ji96  = ε̇(Ji96,  σ_ax.*1e6, 1473., 2e9, fH2O)

# res = 800
# fig = Figure(resolution = (res, res), fontsize=25)
# ax = Axis(fig[1, 1], title = L"$$Garnet creep laws", xlabel = L"$\sigma$ [MPa]", ylabel = L"$\dot{\varepsilon}$ [s$^{-1}$]", aspect=1.0, xscale=log10, yscale=log10)
# lines!(ax, σ_ax, ε̇_Xu13,  label="Xu et al. (2013)")
# lines!(ax, σ_ax, ε̇_Ji96,  label="Ji & Martignole (1996)")
# lines!(ax, σ_ax, ε̇_Mei10, label="Xu et al. (2013)")
# lines!(ax, σ_ax, ε̇_Li06,  label="Li et al. (2006)")
# axislegend("References", position = :rb)
# display(fig)

# if printfig save( string(@__DIR__,  "/Figure8_Xu13.png"), fig, px_per_unit = 4) end

# #--------------------------------------------------------------#
# # Reproduce figure from Zulauf et al 202+
# e      = 1e-13
# T_ax   = LinRange(100, 500, 400) .+ 273.15 
# P_ax   = LinRange(0.0, 1.9, 400) .* 1e9 
# σ_Wu13 = σ(Xu13,  e, T_ax, P_ax', fH2O)

# σ_brit = Coh.*cos(phi) .+ P_ax'.*sin(phi) .+ 0 .* T_ax

# σ_Wu13[σ_Wu13.>σ_brit] .= σ_brit[σ_Wu13.>σ_brit]
# T_Wu13 = T_BDT(Xu13,  e, P_ax, fH2O) .- 273.15

# res = 800
# fig = Figure(resolution = (res, res), fontsize=25)
# ax = Axis(fig[1, 1], title = L"$$Garnet deformation map (Xu et al., 2013, $f_\textrm{H2O}$ = 12 GPa)", xlabel = L"$T$ [°C]", ylabel = L"$P$ [GPa]", aspect=1.0)
# heatmap!(ax, T_ax.-273.15, P_ax./1e9, log10.(σ_Wu13) , colormap=(:lajolla))
# contour!(ax, T_ax.-273.15, P_ax./1e9, σ_Wu13, levels=[800e6, 400e6, 200e6, 100e6], color=:white)
# lines!(ax, T_Wu13, P_ax./1e9, color=:white, linewidth=10 )
# text!(ax, 300, 0.7, text=L"$$Brittle-Ductile transition", rotation=π/2-0.15, color=:white)
# text!(ax, 400, 0.7, text=L"$100$ MPa", rotation=π/2-0.33, color=:white)
# text!(ax, 380, 0.9, text=L"$200$ MPa", rotation=π/2-0.33, color=:white)
# text!(ax, 375, 1.3, text=L"$400$ MPa", rotation=π/2-0.33, color=:white)
# text!(ax, 365, 1.6, text=L"$800$ MPa", rotation=π/2-0.30, color=:white)
# if printfig save( string(@__DIR__,  "/GarnetDeformationMap_Xu13.png"), fig, px_per_unit = 4) end
# display(fig)


# #--------------------------------------------------------------#
# # Reproduce figure from Zulauf et al 202+
# e      = 1e-13
# T_ax   = LinRange(100, 850, 400) .+ 273.15 
# P_ax   = LinRange(0.0, 1.9, 400) .* 1e9 
# σ_Li06 = σ(Li06,  e, T_ax, P_ax', fH2O)

# σ_brit = Coh.*cos(phi) .+ P_ax'.*sin(phi) .+ 0 .* T_ax

# σ_Li06[σ_Li06.>σ_brit] .= σ_brit[σ_Li06.>σ_brit]
# T_Li06 = T_BDT(Li06,  e, P_ax, fH2O) .- 273.15

# res = 800
# fig = Figure(resolution = (res, res), fontsize=25)
# ax = Axis(fig[1, 1], title = L"$$Garnet deformation map (Li et al., 2006)", xlabel = L"$T$ [°C]", ylabel = L"$P$ [GPa]", aspect=1.0)
# heatmap!(ax, T_ax.-273.15, P_ax./1e9, log10.(σ_Li06) , colormap=(:lajolla))
# contour!(ax, T_ax.-273.15, P_ax./1e9, σ_Li06, levels=[800e6, 400e6, 200e6, 100e6], color=:white)
# lines!(ax, T_Li06, P_ax./1e9, color=:white, linewidth=10 )
# # poly!(ax, [Rectangle{Float64}(125, 125, 50, 50)], color = :gray, strokewidth = 2, strokecolor = :red)
# # text!(ax, 300, 0.7, text=L"$$Brittle-Ductile transition", rotation=π/2-0.15, color=:white)
# # text!(ax, 400, 0.7, text=L"$100$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 380, 0.9, text=L"$200$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 375, 1.3, text=L"$400$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 365, 1.6, text=L"$800$ MPa", rotation=π/2-0.30, color=:white)
# if printfig save( string(@__DIR__,  "/GarnetDeformationMap_Li06.png"), fig, px_per_unit = 4) end
# display(fig)

# #--------------------------------------------------------------#
# # Reproduce figure from Zulauf et al 202+
# e      = 1e-13
# T_ax   = LinRange(100, 850, 400) .+ 273.15 
# P_ax   = LinRange(0.0, 1.9, 400) .* 1e9 
# σ_Ji96 = σ(Ji96,  e, T_ax, P_ax', fH2O)

# σ_brit = Coh.*cos(phi) .+ P_ax'.*sin(phi) .+ 0 .* T_ax

# σ_Ji96[σ_Ji96.>σ_brit] .= σ_brit[σ_Ji96.>σ_brit]
# T_Ji96 = T_BDT(Ji96,  e, P_ax, fH2O) .- 273.15

# res = 800
# fig = Figure(resolution = (res, res), fontsize=25)
# ax = Axis(fig[1, 1], title = L"$$Garnet deformation map (Li et al., 2006)", xlabel = L"$T$ [°C]", ylabel = L"$P$ [GPa]", aspect=1.0)
# heatmap!(ax, T_ax.-273.15, P_ax./1e9, log10.(σ_Ji96) , colormap=(:lajolla))
# contour!(ax, T_ax.-273.15, P_ax./1e9, σ_Ji96, levels=[800e6, 400e6, 200e6, 100e6], color=:white)
# lines!(ax, T_Ji96, P_ax./1e9, color=:white, linewidth=10 )
# # poly!(ax, [Rectangle{Float64}(125, 125, 50, 50)], color = :gray, strokewidth = 2, strokecolor = :red)
# # text!(ax, 300, 0.7, text=L"$$Brittle-Ductile transition", rotation=π/2-0.15, color=:white)
# # text!(ax, 400, 0.7, text=L"$100$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 380, 0.9, text=L"$200$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 375, 1.3, text=L"$400$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 365, 1.6, text=L"$800$ MPa", rotation=π/2-0.30, color=:white)
# if printfig save( string(@__DIR__,  "/GarnetDeformationMap_Ji96.png"), fig, px_per_unit = 4) end
# display(fig)

# #--------------------------------------------------------------#
# # Reproduce figure from Zulauf et al 202+
# e      = 1e-13
# T_ax   = LinRange(100, 850, 400) .+ 273.15 
# P_ax   = LinRange(0.0, 1.9, 400) .* 1e9 
# σ_Mei10 = σ(Mei10,  e, T_ax, P_ax', fH2O)

# σ_brit = Coh.*cos(phi) .+ P_ax'.*sin(phi) .+ 0 .* T_ax

# σ_Mei10[σ_Mei10.>σ_brit] .= σ_brit[σ_Mei10.>σ_brit]
# T_Mei10 = T_BDT(Mei10,  e, P_ax, fH2O) .- 273.15

# res = 800
# fig = Figure(resolution = (res, res), fontsize=25)
# ax = Axis(fig[1, 1], title = L"$$Garnet deformation map (Li et al., 2006)", xlabel = L"$T$ [°C]", ylabel = L"$P$ [GPa]", aspect=1.0)
# heatmap!(ax, T_ax.-273.15, P_ax./1e9, log10.(σ_Mei10) , colormap=(:lajolla))
# contour!(ax, T_ax.-273.15, P_ax./1e9, σ_Mei10, levels=[800e6, 400e6, 200e6, 100e6], color=:white)
# lines!(ax, T_Mei10, P_ax./1e9, color=:white, linewidth=10 )
# # poly!(ax, [Rectangle{Float64}(125, 125, 50, 50)], color = :gray, strokewidth = 2, strokecolor = :red)
# # text!(ax, 300, 0.7, text=L"$$Brittle-Ductile transition", rotation=π/2-0.15, color=:white)
# # text!(ax, 400, 0.7, text=L"$100$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 380, 0.9, text=L"$200$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 375, 1.3, text=L"$400$ MPa", rotation=π/2-0.33, color=:white)
# # text!(ax, 365, 1.6, text=L"$800$ MPa", rotation=π/2-0.30, color=:white)
# if printfig save( string(@__DIR__,  "/GarnetDeformationMap_Mei10.png"), fig, px_per_unit = 4) end
# display(fig)

end 

main()
