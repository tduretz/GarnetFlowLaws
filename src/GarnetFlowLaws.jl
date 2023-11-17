module GarnetFlowLaws

# Local dependency
using Enzyme

# Constants
const R   = 8.316 # Gas constant
export R

# Anonymous functions
T_BDT(p, Eii, P, fH2O, λf, phi, Coh) =  (p.Q .+ P .* p.V) ./ (R .* log.(p.A .* fH2O .^ p.r .* (2*(Coh .* cos(phi) .+ P .* sin(phi)*λf)) .^ p.n ./ Eii))

ε̇(p, σ, T, P, fH2O)  = p.A*σ.^p.n*fH2O.^p.r*exp(- (p.Q + P*p.V)/R/T)

σ(p, ε̇, T, P, fH2O) = (ε̇ .* fH2O .^ (-p.r) .* exp.((p.Q .+ P .* p.V) ./ (R .* T)) ./ p.A) .^ (1 ./ p.n)

export T_BDT, ε̇, σ

# Regular functions
function T_BDT2(p, Eii, P, λf, phi, Coh) 
    fH2O   = 1e8
    T      = 1000
    for iter=1:500
        T = T_BDT(p, Eii, P, fH2O, λf, phi, Coh)
        fH2O =  PSfugacity(P, T)
    end
    return T
end

function Compute_c!( c, ci, T)
    for i in axes(ci,1)
        c[i]  =  ci[i,1]*T^(-4) + ci[i,2]*T^(-2) + ci[i,3]*T^(-1) + ci[i,4] + ci[i,5]*T + ci[i,6]*T^2
    end
end

function PSeos(V, T, targetP, c)  # cc/mol, Kelvins, Pa
    R  =  8314510  # Pa.cc/K/mol
    ρ  =  1/V  # mol/cc
    P = (ρ + c[1]*ρ^2-ρ^2*((c[3]+2*c[4]*ρ+3*c[5]*ρ^2
              +4*c[6]*ρ^3)/(c[2]+c[3]*ρ+c[4]*ρ^2+c[5]*ρ^3
              +c[6]*ρ^4)^2)+c[7]*ρ^2*exp(-c[8]*ρ)
              +c[9]*ρ^2*exp(-c[10]*ρ))*R*T
    return P-targetP  # Pa
end
      
function PSV(P, T, c)  # Pa, Kelvins
    V  =  10.
    f0 = 0.
    for iter = 1:100
        f    =  PSeos(V, T, P, c)
        ∂f∂v =  autodiff(Enzyme.Reverse, PSeos, Active, Active(V), Const(T), Const(P), Const(c))[1][1]
        if iter==1 f0 = f end
        if abs(f)<1e-10 || abs(f)/abs(f0)<1e-10 break end
        # @show (iter, f, f/f0, V)
        V   -=  f/∂f∂v
    end
    return V
end
   
function PSfugacity(P, T)  # Pa, Kelvins
    # H20
    ci  =  [
        0.          0.         0.246576e6 0.513599e+2  0.           0.
        0.          0.         0.586389e0 -.286469e-2  0.313755e-4  0.
        0.          0.         -.627838e1 0.147959e-1  0.357795e-3  0.154329e-7
        0.          0.         0.         -.427198e+0  -.163251e-4  0.
        0.          0.         0.566549e4 -.165801e+2  0.765607e-1  0.
        0.          0.         0.         0.109117e+0  0.           0.
        0.388786e13 -.134948e9 0.309165e6 0.755911e+1  0.           0.
        0.          0.         -.655378e5 0.188106e+3  0.           0.
        -.141824e14 0.181653e9 -.197690e6 -.235303e+2  0.           0.
        0.          0.         0.920993e5 0.122467e+3  0.           0.
    ]
    R = 8314510  # Pa.cc/K/mol
    c = zeros(10)
    Compute_c!( c, ci, T)
    V = PSV(P, T, c)
    ρ = 1/V  # mol/cc
    if ρ<0. ρ=1e-10 end 
    fug = exp(log(ρ) + c[1]*ρ + (1/(c[2]+c[3]*ρ+c[4]*ρ^2
                                         +c[5]*ρ^3+c[6]*ρ^4)-1/c[2])
                 -c[7]/c[8]*(exp(-c[8]*ρ)-1)
                 -c[9]/c[10]*(exp(-c[10]*ρ)-1)
                 +P/(ρ*R*T)
                 +log(R*T)-1)
    return fug  
end

export T_BDT2, PSfugacity, PSV, Compute_c!

# Data
include("GarnetCreepLaws.jl")
export Ji96, Mei10, Xu13, Li06, Ka95

end # module GarnetFlowLaws
