using Pkg

# Pkg.add("DifferentialEquations")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Distributions")
# Pkg.add("LinearAlgebra")
# Pkg.add("Tables")
# Pkg.add("Random")

using DifferentialEquations
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Tables
using Dates
using Random

idx = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])

dist_df = DataFrame(CSV.File("files/coord_distances_R3.csv"))
ll_df = DataFrame(CSV.File("files/morphotype_dist_data.csv"))
ll_df[:,:n_trees] = convert.(Int,round.(ll_df.n_trees, digits = 0))

param_df = DataFrame(CSV.File("files/nt4s_top30.csv"))

model  = "no trees"
dir = "realization_op/op_nt1"
get_all_data = false

a = now()
println("Time started: "*string(a)*" hours")

const d = 0.2
const eps_val = 0.1
const gen = 200
const index = 30

const bi = 1
const br = 1

const r = 0.2
const gamma = 0.2
const μ1 = 0.55
const μ2 = 0.55
const δ1 = 0.083
const δ2 = 0.083
const k1 = Int(20)
const k2 = Int(20)

const νSDO = 0.1
const νSGR = 0.1
const νMDO = 0.1
const νMGR = 0.1

if model == "no trees"
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

global number_of_pops = length(Tables.matrix(unique(select(dist_df, "num"))))::Int
#number_of_pops = 1

global pop_nums = (1:number_of_pops)

function twostrain_SEIR(du,u,p,t)
    μ1, μ2, δ1, δ2, k1, k2, C1, C2, rho = p

    du[1] = -(u[4 + k1 + k2]*u[2]*u[1]) - (u[5 + k1 + k2]*u[3]*u[1])                     # S
    du[2] = -(C1^2*u[2]^2*u[4 + k1 + k2]) - (u[5 + k1 + k2]*rho*C1*C2*u[2]*u[3])            # nu_bar 1
    du[3] = -(C2^2*u[3]^2*u[5 + k1 + k2]) - (u[4 + k1 + k2]*rho*C1*C2*u[2]*u[3])            # nu_bar2

    du[4] = u[4 + k1 + k2]*u[2]*u[1] - k1*δ1*u[4]                                        # E1 first

    for i in 2:k1
        du[4 + i - 1] = (k1*δ1*u[4 + i - 2]) - k1*δ1*u[4 + i - 1]                      # E1 second
    end

    du[4 + k1] = u[5 + k1 + k2]*u[3]*u[1] - k2*δ2*u[4 + k1]                              # E2 first

    for i in 2:k2
        du[4 + k1 + i - 1] = k2*δ2*u[4 + k1 + i - 2] - k2*δ2*u[4 + k1 + i - 1]       # E2 second
    end

    du[4 + k1 + k2] = k1*δ1*u[4 + k1 - 1] - μ1*u[4 + k1 + k2]                        # P1
    du[5 + k1 + k2] = k2*δ2*u[4 + k1 + k2 - 1] - μ2*u[5 + k1 + k2]                   # P2
    du[6 + k1 + k2] = k1*δ1*u[4 + k1 - 1]                                            # psum1
    du[7 + k1 + k2] = k2*δ2*u[4 + k1 + k2 - 1]                                       # psum2
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

    frac_i1 = (I1_end)/(S_old)
    frac_i2 = (I2_end)/(S_old)
    Nt1 = stoch*S_end*(r + r*(si*νi + sr*νr))

    n1t1 = (1/(1 + si*νi + sr*νr))*(νi + si*(νi^2)*((bi^2)*(Ci^2) + 1) + sr*νr*νi*(ρ*bi*br*Ci*Cr + 1))
    n2t1 = (1/(1 + si*νi + sr*νr))*(νr + sr*(νr^2)*((br^2)*(Cr^2) + 1) + si*νr*νi*(ρ*bi*br*Ci*Cr + 1))

    Z1t1 = ϕ1*I1_end + gamma*Z1old
    Z2t1 = ϕ2*I2_end + gamma*Z2old

    ow_list = [frac_i1, frac_i2, Nt1, Z1t1, Z2t1, n1t1, n2t1]

    return ow_list
end

function place_tree_sp(n_do, n)
    if n_do >= 1
        tree_nums = sample(1:n,n_do, replace = false)
        tree_arr = repeat(["GR"],outer = 37)
        tree_arr[tree_nums] .= "DO"
    else
        tree_arr = repeat(["GR"],outer = 37)
    end
    return tree_arr
end

function get_prop(ϵ, distance_matrix)
    exp_mat = exp.(-ϵ*distance_matrix)
    exp_mat[diagind(exp_mat)] .= 0.00
    exp_row_sum = sum(exp_mat,dims = 2)
    prop_mat = exp_mat./exp_row_sum
    return(prop_mat)
end

function move_pops(pop_mat, t, prop_disp, d)
    B = Array{Float64}(undef, 1, number_of_pops)
    for i in 1:number_of_pops
        B[1,i] = pop_mat[i,t]*(1 - d) + sum(d.*pop_mat[:,t].*prop_disp[:,i])
    end

    return(B)
end

function move_pops_nu(S_mat, nu_mat, t, prop_disp, d)
    B = Array{Float64}(undef, 1, number_of_pops)
    for i in 1:number_of_pops
        B[1,i] = (S_mat[i,t]*nu_mat[i,t]*(1 - d) + sum(d.*S_mat[:,t].*prop_disp[:,i].*nu_mat[:,t]))/(S_mat[i,t]*(1-d) + sum(d.*S_mat[:,t].*prop_disp[:,i]))
    end
    return(B)
end

dist_mat = Array{Float64}(undef, number_of_pops, number_of_pops)
for i in pop_nums
    temp = select(filter(:ref=>==(i),dist_df), [:distance])
    dist_mat[:,i] = transpose(convert(Vector{Float64}, temp.distance))
end

function dbinom(n,k,p)
    b = Binomial(k,p)
    return pdf(b,n)
end

function run_simulation(pdoug::Int64, phiS::Float64, phiM::Float64, rep::Int64, tree_vals, sS::Float64,sM::Float64,rho::Float64,sigma::Float64)

    if number_of_pops >2
        prop_disp = get_prop(eps_val,dist_mat)::Matrix{Float64}
    elseif number_of_pops == 2
        prop_disp = [0 1; 1 0]::Matrix{Float64}
    else
        prop_disp = Array{Float64}(undef, 1, number_of_pops)
        prop_disp[1,1] = 0
    end

    S_pop = zeros(number_of_pops, gen+1)::Matrix{Float64}
    Z1_pop = zeros(number_of_pops, gen+1)::Matrix{Float64}
    Z2_pop = zeros(number_of_pops, gen+1)::Matrix{Float64}
    frac_pop1 = zeros(number_of_pops, gen+1)::Matrix{Float64}
    frac_pop2 = zeros(number_of_pops, gen+1)::Matrix{Float64}
    nu1_track = zeros(number_of_pops, gen+1)::Matrix{Float64}
    nu2_track = zeros(number_of_pops, gen+1)::Matrix{Float64}

    initS = rand(Uniform(0.4,1.8),1)[1]::Float64
    initMNPV = rand(LogNormal(-1,0.5),1)[1]::Float64
    initSNPV = rand(LogNormal(-1,0.5),1)[1]::Float64

    for p in 1:number_of_pops

        S_pop[p,1] = initS
        Z1_pop[p,1] = initSNPV
        Z2_pop[p,1] = initMNPV

        if tree_vals[p] == "DO"
            ν_SNPV = νSDO # ν_SNPV_DO
            ν_MNPV = νMDO
        else
            ν_SNPV = νSGR
            ν_MNPV = νMGR
        end

        nu1_track[p,1] = ν_SNPV
        nu2_track[p,1] = ν_MNPV

    end

    for t in 2:(gen + 1)

        S_move = move_pops(S_pop,t-1, prop_disp, d)
        Z1_move = move_pops(Z1_pop,t-1, prop_disp, d)
        Z2_move = move_pops(Z2_pop,t-1, prop_disp, d)
        nu1_move = move_pops_nu(S_pop, nu1_track,t-1, prop_disp, d)
        nu2_move = move_pops_nu(S_pop, nu2_track,t-1, prop_disp, d)

        en = rand(Normal(0,sigma),1)
        stoch = exp(en[1])

        for p in 1:number_of_pops

            if tree_vals[p] == "DO"
                C_SNPV = C_SNPV_DO
                C_MNPV = C_MNPV_DO

            else
                C_SNPV = C_SNPV_GR
                C_MNPV = C_MNPV_GR
            end

            S_init = S_move[p]::Float64
            Z1_init = Z1_move[p]::Float64 #+ 1e-8
            Z2_init = Z2_move[p]::Float64 #+ 1e-8
            nu1_init = nu1_move[p]::Float64
            nu2_init = nu2_move[p]::Float64

            output = run_ode(S_init,Z1_init,Z2_init,nu1_init,nu2_init,C_SNPV,C_MNPV,rho)

            simulation = output[1]
            num_equations = output[2]
            
            end_t = length(simulation.t)::Int64

            S_end = simulation[1,end_t]::Float64
            nu1_end = simulation[2,end_t]::Float64
            nu2_end = simulation[3,end_t]::Float64
            I1_end = simulation[num_equations - 1,end_t]::Float64
            I2_end = simulation[num_equations,end_t]::Float64

            #S_end, I1_end, I2_end, νi, νr, Ci, Cr, bi, br, sr, si, r, ϕ1, ϕ2, γ, ρ
            over_winter = inter_annual_evolution(S_end, I1_end, I2_end, nu1_end, nu2_end, S_init, Z1_init, Z2_init, C_SNPV, C_MNPV, sS, sM, phiS, phiM, rho, stoch)::Vector{Float64}

            frac_pop1[p,t] = over_winter[1]
            frac_pop2[p,t] = over_winter[2]
            S_pop[p,t] = over_winter[3]
            Z1_pop[p,t] = over_winter[4]
            Z2_pop[p,t] = over_winter[5]
            nu1_track[p,t] = over_winter[6]
            nu2_track[p,t] = over_winter[7]

        end
    end

    if get_all_data == true

        S_means = mean(S_pop,dims = 1)
        Z1_means = mean(Z1_pop,dims = 1)
        Z2_means = mean(Z2_pop,dims = 1)
        nu1_means = mean(nu1_track,dims = 1)
        nu2_means = mean(nu2_track,dims = 1)
        fraci1_means = mean(frac_pop1,dims = 1)
        fraci2_means = mean(frac_pop2,dims = 1)

        S_dat = DataFrame(Tables.table(transpose(S_means)))

        rows = size(S_dat)[1]
        ncols = size(S_dat)[2]

        insertcols!(S_dat, 1, :Year => 1:rows)
        insertcols!(S_dat, 1, :Species => "S")

        Z1_dat = DataFrame(Tables.table(transpose(Z1_means)))
        insertcols!(Z1_dat, 1, :Year => 1:rows)
        insertcols!(Z1_dat, 1, :Species => "Z1")

        Z2_dat = DataFrame(Tables.table(transpose(Z2_means)))
        insertcols!(Z2_dat, 1, :Year => 1:rows)
        insertcols!(Z2_dat, 1, :Species => "Z2")

        fracI1_dat = DataFrame(Tables.table(transpose(fraci1_means)))
        insertcols!(fracI1_dat, 1, :Year => 1:rows)
        insertcols!(fracI1_dat, 1, :Species => "FracI1")

        fracI2_dat = DataFrame(Tables.table(transpose(fraci2_means)))
        insertcols!(fracI2_dat, 1, :Year => 1:rows)
        insertcols!(fracI2_dat, 1, :Species => "FracI2")

        nu1_dat = DataFrame(Tables.table(transpose(nu1_means)))
        insertcols!(nu1_dat, 1, :Year => 1:rows)
        insertcols!(nu1_dat, 1, :Species => "nu1")

        nu2_dat = DataFrame(Tables.table(transpose(nu2_means)))
        insertcols!(nu2_dat, 1, :Year => 1:rows)
        insertcols!(nu2_dat, 1, :Species => "nu2")

        all_long = vcat(S_dat, Z1_dat, Z2_dat, nu1_dat, nu2_dat,fracI1_dat, fracI2_dat)

        all_long.rep .= rep
        all_long.pdoug .= pdoug
        all_long.sS_DO .= sS_DO
        all_long.sS_GR .= sS_GR
        all_long.sM_DO .= sM_DO
        all_long.sM_GR .= sM_GR
        all_long.phi .= phiS
        all_long.rho .= rho
        all_long.sigma .= sigma
    
    end

    total_frac = (frac_pop1 + frac_pop2)::Matrix{Float64}
    total_frac_means = mean(total_frac,dims = 1)::Matrix{Float64}

    frac2_means = mean(frac_pop2,dims = 1)::Matrix{Float64}

    mnpv_frac = (frac2_means ./ total_frac_means)::Matrix{Float64}

    f_index = findall(x -> x >= index/100, total_frac_means[1,:])::Vector{Int64}

    m_frac_mean = mean(mnpv_frac[:,f_index])::Float64

    mean_last_S = mean(S_pop[:,(gen-50):gen])::Float64
    max_last_S = maximum(S_pop[:,(gen-50):gen])::Float64
    min_last_S = minimum(S_pop[:,(gen-50):gen])::Float64

    mean_last_Z1 = mean(Z1_pop[:,(gen-50):gen])::Float64
    mean_last_Z2 = mean(Z2_pop[:,(gen-50):gen])::Float64

    if mean_last_S <= 12 && mean_last_S >= 1e-2 && max_last_S <= 200 && min_last_S >= 1e-10
        extinct_S = 0::Int
    else
        extinct_S = 1::Int
    end

    if mean_last_Z1 <= 1e-10 && extinct_S == 0
        extinct_SNPV = 1::Int  
    else
        extinct_SNPV = 0::Int
    end

    if mean_last_Z2 <= 1e-10 && extinct_S == 0
        extinct_MNPV = 1::Int
    else
        extinct_MNPV = 0::Int
    end

    if extinct_S == 1 #|| extinct_SNPV == 1 || extinct_MNPV == 1
        qual = 0::Int
    else
        qual = 1::Int
    end

    small_df = filter(:n_trees => ==(pdoug), ll_df)::DataFrame

    probs_df = DataFrame()
    for i in 1:length(small_df.id)
        if qual == 1
            if isnan(m_frac_mean) == true
                p1 = NaN
            else
                p1 = dbinom(small_df.MNPV[i],small_df.total[i],m_frac_mean)::Float64
            end
        else
            p1 = 0 # 1e-150::Float64
        end
 
        temp_df = DataFrame(LL = p1, id = small_df.id[i])

        probs_df = vcat(probs_df,temp_df)
    end

    probs_df.quality .= qual
    probs_df.rep .= rep
    probs_df.pdoug .= pdoug

    probs_df.extinct_S .= extinct_SNPV
    probs_df.extinct_M .= extinct_MNPV
    probs_df.mean_pMNPV .= m_frac_mean
    probs_df.sS_DO .= sS_DO
    probs_df.sS_GR .= sS_GR
    probs_df.sM_DO .= sM_DO
    probs_df.sM_GR .= sM_GR
    probs_df.phi .= phiS
    probs_df.rho .= rho
    probs_df.sigma .= sigma

    if get_all_data == false
        return probs_df 
    else 
        return probs_df, all_long
    end

end

global sS = [param_df.sS[idx]]
global sM = [param_df.sM[idx]]
global phi_S = global phi_M = [param_df.phi[idx]]
global rho = [param_df.rho[idx]]

#length(sr_test)*length(si_test)*200*1*3.9/(60*60)*1.2

global ndoug = 18
global sigma = 1.25

if isdir(dir) == false
    mkdir(dir) 
end

if isdir(dir*"ll/") == false
    if get_all_data == false
        mkdir(dir*"ll/")
    else
        mkdir(dir*"ll/")
        mkdir(dir*"all/")
    end
end

ll_data = DataFrame()
#all_data = DataFrame()

for p in 1:length(phi_S)
    for rp in 1:10 # change to 2000 if not testing

        tree_vals = place_tree_sp(ndoug,number_of_pops)
        
        rep = rp
        
        output_sim = run_simulation(ndoug,phi_S[p],phi_M[p],rep,tree_vals,sS[p],sM[p],rho[p],sigma)

        if get_all_data == false
            ll_data = vcat(ll_data,output_sim)
        else
            ll_data = vcat(ll_data,output_sim[1])
            all_data = vcat(all_data,output_sim[2])
        end

    end
end

sim_info = "pset_"*string(idx)*"_test.csv"

if get_all_data == false
    CSV.write(dir*"ll/"*sim_info, ll_data, header = true)
else
    CSV.write(dir*"ll/"*sim_info, ll_data, header = true)
    CSV.write(dir*"all/"*sim_info, all_data, header = true)
end

b = now()
time_elapsed = round((b-a).value/(60000*60),digits = 2)
println("Time elapsed: "*string(time_elapsed)*" hours")