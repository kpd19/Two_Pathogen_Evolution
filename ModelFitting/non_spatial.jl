using DifferentialEquations
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Tables
using Dates
using Random

# running array jobs on the midway computing cluster, idx is the index of the array
trees = ["DO","GR"]

model = "host trees"
dir = "realization_op/non_spatial_evolution/"

a = now()
println("Time started: "*string(a)*" hours")

const gen = 200 # number of generations
const index = 30 # 0.3 index, only

const bi = 1 # heritability for pathogen 1
const br = 1 # heritability for pathogen 2

const r = 0.2 # baseline host reproduction rate
const gamma = 0.2 # between season pathogen 

const μ1 = 0.55 # decay rate of SNPV
const μ2 = 0.55 # decay rate of MNPV 
const δ1 = 0.083 # rate of movement between exposed classes for SNPV
const δ2 = 0.083 # rate of movement between exposed classes for MNPV
const k1 = Int(20) # number of exposed classes for SNPV
const k2 = Int(20) # number of exposed classes for MNPV 

const νSDO = 0.1 # initial value for transmission rate for SNPV on Douglas-fir
const νSGR = 0.1 # SNPV on grand fir
const νMDO = 0.1 # MNPV on Douglas-fir
const νMGR = 0.1 # MNPV on grand fir

if model == "host trees"
    const C_SNPV_DO = 3.7 # Heterogeneity in infectiousness for SNPV on Douglas-fir
    const C_SNPV_GR = 3.7 # Heterogeneity in infectiousness for SNPV on grand fir
    const C_MNPV_DO = 2.5 # MNPV on Douglas-fir
    const C_MNPV_GR = 3.3 # MNPV on grand fir
elseif model == "no trees"
    const C_SNPV_DO = 3.4 # model 3 no tree species
    const C_SNPV_GR = 3.4 # model 3 no tree species
    const C_MNPV_DO = 2.9 # model 3 no tree species
    const C_MNPV_GR = 2.9 # model 3 no tree species
elseif model == "no morphotype"
    const C_SNPV_DO = 3.1 # model 3 no tree species or morph
    const C_SNPV_GR = 3.1 # model 3 no tree species or morph
    const C_MNPV_DO = 3.1 # model 3 no tree species or morph
    const C_MNPV_GR = 3.1 # model 3 no tree species or morph
else 
    print("Error: No model specified")
end

function twostrain_SEIR(du,u,p,t)
    μ1, μ2, δ1, δ2, k1, k2, C1, C2, rho = p

    du[1] = -(u[4 + k1 + k2]*u[2]*u[1]) - (u[5 + k1 + k2]*u[3]*u[1])                        # S
    du[2] = -(C1^2*u[2]^2*u[4 + k1 + k2]) - (u[5 + k1 + k2]*rho*C1*C2*u[2]*u[3])            # nu_bar 1
    du[3] = -(C2^2*u[3]^2*u[5 + k1 + k2]) - (u[4 + k1 + k2]*rho*C1*C2*u[2]*u[3])            # nu_bar2

    du[4] = u[4 + k1 + k2]*u[2]*u[1] - k1*δ1*u[4]                                           # E1 first

    for i in 2:k1
        du[4 + i - 1] = (k1*δ1*u[4 + i - 2]) - k1*δ1*u[4 + i - 1]                           # E1 second
    end

    du[4 + k1] = u[5 + k1 + k2]*u[3]*u[1] - k2*δ2*u[4 + k1]                                 # E2 first

    for i in 2:k2
        du[4 + k1 + i - 1] = k2*δ2*u[4 + k1 + i - 2] - k2*δ2*u[4 + k1 + i - 1]              # E2 second
    end

    du[4 + k1 + k2] = k1*δ1*u[4 + k1 - 1] - μ1*u[4 + k1 + k2]                               # P1
    du[5 + k1 + k2] = k2*δ2*u[4 + k1 + k2 - 1] - μ2*u[5 + k1 + k2]                          # P2
    du[6 + k1 + k2] = k1*δ1*u[4 + k1 - 1]                                                   # psum1
    du[7 + k1 + k2] = k2*δ2*u[4 + k1 + k2 - 1]                                              # psum2
    nothing
end

function run_ode(S0, p0_1, p0_2, ν1, ν2, C1, C2, rho)
    rho = rho

    e1 = zeros(k1)
    e2 = zeros(k2)

    u0 = hcat(S0, ν1, ν2,transpose(e1),transpose(e2),p0_1, p0_2, 0,0)
    p = (μ1, μ2, δ1, δ2, k1, k2, C1, C2, rho)

    num_eq = length(u0)

    tspan = (0.0,70.0)
    prob = ODEProblem{true}(twostrain_SEIR,u0,tspan,p)
    sol = solve(prob, Tsit5(), abstol = 1e-8, reltol = 1e-8, isoutofdomain=(u,p,t) -> any(x -> x < 0, u), save_everystep = false, verbose = false)

    return(sol,num_eq)

end

function inter_annual_evolution(S_end, I1_end, I2_end, νi, νr, S_old, Z1old, Z2old, Ci, Cr, si, sr, ϕ1, ϕ2, ρ, stoch)

    frac_i1 = (I1_end)/(S_old) # fraction infected by SNPV
    frac_i2 = (I2_end)/(S_old) # fraction infected by MNPV
    Nt1 = stoch*S_end*(r + r*(si*νi + sr*νr)) # Host population after overwintering

    n1t1 = (1/(1 + si*νi + sr*νr))*(νi + si*(νi^2)*((bi^2)*(Ci^2) + 1) + sr*νr*νi*(ρ*bi*br*Ci*Cr + 1)) # transmission risk after overwintering for SNPV
    n2t1 = (1/(1 + si*νi + sr*νr))*(νr + sr*(νr^2)*((br^2)*(Cr^2) + 1) + si*νr*νi*(ρ*bi*br*Ci*Cr + 1)) # transmission risk after overwintering for MNPV

    Z1t1 = ϕ1*I1_end + gamma*Z1old # SNPV populationa after overwintering
    Z2t1 = ϕ2*I2_end + gamma*Z2old # MNPV population after overwintering

    ow_list = [frac_i1, frac_i2, Nt1, Z1t1, Z2t1, n1t1, n2t1]

    return ow_list
end


function run_simulation(tree_sp, phiS::Float64, phiM::Float64, rep::Int64, sS_DO::Float64,sS_GR::Float64,sM_DO::Float64,sM_GR::Float64,rho::Float64,sigma::Float64)

    S_pop = zeros(gen+1)::Vector{Float64} # host
    Z1_pop = zeros(gen+1)::Vector{Float64} # SNPV
    Z2_pop = zeros(gen+1)::Vector{Float64} # MNPV
    frac_pop1 = zeros(gen+1)::Vector{Float64} # fraction infected SNPV
    frac_pop2 = zeros(gen+1)::Vector{Float64} # fraction infected MNPV
    nu1_track = zeros(gen+1)::Vector{Float64} # transmission risk SNPV
    nu2_track = zeros(gen+1)::Vector{Float64} # transmission risk MNPV

    initS = rand(Uniform(0.4,1.8),1)[1]::Float64 # random initial conditions for host
    initMNPV = rand(LogNormal(-1,0.5),1)[1]::Float64 # random initial conditions for SNPV
    initSNPV = rand(LogNormal(-1,0.5),1)[1]::Float64 # random initial conditions for MNPV

    S_pop[1] = initS # first generation host
    Z1_pop[1] = initSNPV # first generation SNPV
    Z2_pop[1] = initMNPV # first generation MNPV

    if tree_sp == "DO"
        ν_SNPV = νSDO # transmission risk SNPV on Douglas-fir
        ν_MNPV = νMDO # tranmission risk MNPV on Douglas-fir
    else
        ν_SNPV = νSGR # transmission risk SNPV on grand fir
        ν_MNPV = νMGR # tranmission risk MNPV on grand fir
    end

    nu1_track[1] = ν_SNPV # first generation SNPV transmission risk
    nu2_track[1] = ν_MNPV # first generation MNPV transmission risk

    for t in 2:(gen + 1)

        en = rand(Normal(0,sigma),1) # environmental stochasticity
        stoch = exp(en[1]) # environmental stochasticity

        if tree_sp == "DO" # setting heterogeneity in transmission risk (C) and tradeoff parameter based on tree species
            C_SNPV = C_SNPV_DO
            C_MNPV = C_MNPV_DO

            sS = sS_DO
            sM = sM_DO
        else
            C_SNPV = C_SNPV_GR
            C_MNPV = C_MNPV_GR

            sS = sS_GR
            sM = sM_GR
        end

        S_init = S_pop[t-1]::Float64 
        Z1_init = Z1_pop[t-1]::Float64 #+ 1e-8
        Z2_init = Z2_pop[t-1]::Float64 #+ 1e-8
        nu1_init = nu1_track[t-1]::Float64
        nu2_init = nu2_track[t-1]::Float64

        output = run_ode(S_init,Z1_init,Z2_init,nu1_init,nu2_init,C_SNPV,C_MNPV,rho) # running ODEs for that generation

        simulation = output[1] # result of ODEs
        num_equations = output[2] # number of equations in ODE set
        
        end_t = length(simulation.t)::Int64 # end of season

        S_end = simulation[1,end_t]::Float64 # host population at end of season
        nu1_end = simulation[2,end_t]::Float64 # tranmission rate of pathogen 1 at end of season
        nu2_end = simulation[3,end_t]::Float64 # transmission rate of pathogen 2 at end of season
        I1_end = simulation[num_equations - 1,end_t]::Float64 # cumulative infected by pathogen 1
        I2_end = simulation[num_equations,end_t]::Float64 # cumulative infected by pathogen 2

        over_winter = inter_annual_evolution(S_end, I1_end, I2_end, nu1_end, nu2_end, S_init, Z1_init, Z2_init, C_SNPV, C_MNPV, sS, sM, phiS, phiM, rho, stoch)::Vector{Float64} # calculate overwintering

        frac_pop1[t] = over_winter[1] # adding new generation to time series 
        frac_pop2[t] = over_winter[2]
        S_pop[t] = over_winter[3]
        Z1_pop[t] = over_winter[4]
        Z2_pop[t] = over_winter[5]
        nu1_track[t] = over_winter[6]
        nu2_track[t] = over_winter[7]

    end

    all_long = DataFrame(time = 1:201, S = S_pop, SNPV = Z1_pop, MNPV = Z2_pop, frac_SNPV = frac_pop1, frac_MNPV = frac_pop2, nu1 = nu1_track, nu2 = nu2_track)

    all_long.rep .= rep
    all_long.tree_sp .= tree_sp
    all_long.sS_DO .= sS_DO
    all_long.sS_GR .= sS_GR
    all_long.sM_DO .= sM_DO
    all_long.sM_GR .= sM_GR
    all_long.phi .= phiS
    all_long.rho .= rho
    all_long.sigma .= sigma
    

    return all_long 
end

global sSD = 1 # tradeoff  parameter for SNPV on Douglas-fir
global sSG = 60 # tradeoff parameter for SNPV on grand fir

global sMD = 6 # tradeoff parameter for MNPV on Douglas-fir
global sMG = 100 # tradeoff parameter for MNPV on grand fir

global phi_S = global phi_M = 25.5 # pathogen overwintering rate
global rho = 0.1 # correlation coefficient between host susceptibilities

global sigma = 0.5 # param_df.sigma #[1e-6, 1e-2, 0.1, 0.5]

if isdir(dir) == false # make directory if it doesn't exitst
    mkdir(dir) 
end

run_simulation("DO",phi_S,phi_S,1,Float64(sSD),Float64(sSG),Float64(sMD),Float64(sMG),rho,sigma)

all_data = DataFrame()

for tree in trees

    for rp in 1:10 # number of realizations
            
        output_sim = run_simulation(tree,Float64(phi_S),Float64(phi_S),rp,Float64(sSD),Float64(sSG),Float64(sMD),Float64(sMG),rho,sigma)

        all_data = vcat(all_data,output_sim)
    end 
end

sim_info = "non_spatial_sig"*string(sigma)*".csv"

CSV.write(dir*sim_info, all_data, header = true)   

ll_data

b = now()
time_elapsed = round((b-a).value/(60000*60),digits = 2)
println("Time elapsed: "*string(time_elapsed)*" hours")