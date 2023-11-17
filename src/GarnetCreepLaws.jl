G  = 8e10      # Shear modulus
Tm = 1853/0.93 # Melting T
    
# Creep - Ji & Matignole (1996)
Ji96 = (
    n  = 2.22,
    Q  = 485e3,
    A  = 5.84e6 *10^(-6*2.22), # from MPa^(-n)/s to Pa^(-n)/s ---------> 2.7952e-07 --- Gut!
    r  = 0.0,
    V  = 0.0
)

# Creep - Mei et al. (2010)
Mei10 = (
    n  = 3.0,
    Q  = 280e3,
    A  = 2.5e8 *10^(-9*3.0), # from GPa^(-n)/s to Pa^(-n)/s
    r  = 0.0,
    V  = 0.0
)

# Creep - Xu et al. (2013)
Xu13 = (
    n  = 3.5,
    r  = 1.0,
    Q  = 215e3,
    A  = 10^(-5.3) *10^(-6*(3.5 + 1.0)), # from MPa^(-n)/s to Pa^(-n)/s
    V  = 28e-6
)

# Creep pyrope - Li et al. 2006
Li06 = (
    n  = 3.2,
    Q  = 270e3,
    A  = exp(15.07) *10^(-9*3.2), # from GPa^(-n)/s to Pa^(-n)/s
    r  = 0.0,
    V  = 0.0
)

# Creep pyrope - Li et al. 2006
Ka95 = (
    n  = 2.7,
    Q  = R*46*Tm,   # g = 46
    A  = exp(46)/(G^2.7),
    r  = 0.0,
    V  = 0.0
)