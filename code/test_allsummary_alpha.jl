
using Pkg
#
# Pkg.add("DifferentialEquations")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Distributions")
# Pkg.add("LinearAlgebra")
# Pkg.add("Tables")

cd("/Users/katherinedixon/Documents/StuffINeed/_Research/Julia_spatial")

using DifferentialEquations
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using Tables

idx = 2
parameters = [10,20,30,40,50,60,70,80]
myparam = parameters[idx + 1]

dist_df = DataFrame(CSV.File("/Users/katherinedixon/Documents/StuffINeed/_Research/_DFTM_2021/Field_2021/_data/coord_distances_R3.csv"))
ll_df = DataFrame(CSV.File("/Users/katherinedixon/Documents/StuffINeed/_Research/Julia_spatial/_data/data_for_ll.csv"))

function twostrain_SEIR(du,u,p,t)
    μ1, μ2, δ1, δ2, k1, k2, C1, C2, alpha = p

    du[1] = -(u[4 + k1 + k2]*u[1]*(u[2] + alpha*u[3])) - (u[5 + k1 + k2]*u[1]*(u[3] + alpha*u[2]))                     # S

    du[2] = -(C1^2*u[2]^2)*(u[4 + k1 + k2] + alpha*u[5 + k1 + k2])       # nu_bar 1
    du[3] = -(C2^2*u[3]^2)*(u[5 + k1 + k2] + alpha*u[4 + k1 + k2])       # nu_bar2

    du[4] = u[4 + k1 + k2]*u[1]*(u[2] + alpha*u[3]) - k1*δ1*u[4]                                        # E1 first

    for i in 2:k1
        du[4 + i - 1] = (k1*δ1*u[4 + i - 2]) - k1*δ1*u[4 + i - 1]                      # E1 second
    end

    du[4 + k1] = u[5 + k1 + k2]*u[1]*(u[3] + alpha*u[2]) - k2*δ2*u[4 + k1]                              # E2 first

    for i in 2:k2
        du[4 + k1 + i - 1] = k2*δ1*u[4 + k1 + i - 2] - k2*δ1*u[4 + k1 + i - 1]       # E2 second
    end

    du[4 + k1 + k2] = k1*δ1*u[4 + k1 - 1] - μ1*u[4 + k1 + k2]                        # P1
    du[5 + k1 + k2] = k2*δ2*u[4 + k1 + k2 - 1] - μ2*u[5 + k1 + k2]                   # P2
    du[6 + k1 + k2] = k1*δ1*u[4 + k1 - 1]                                            # psum1
    du[7 + k1 + k2] = k2*δ2*u[4 + k1 + k2 - 1]                                       # psum2
end

function run_ode(S0, p0_1, p0_2, ν1, ν2, C1, C2, α)
    μ1 = 0.55
    μ2 = 0.55
    δ1 = 0.083
    δ2 = 0.083
    k1 = Int(20)
    k2 = Int(20)
    alpha = α

    e1 = zeros(k1)
    e2 = zeros(k2)

    u0 = hcat(S0, ν1, ν2,transpose(e1),transpose(e2),p0_1, p0_2, 0,0)
    p = (μ1, μ2, δ1, δ2, k1, k2, C1, C2, alpha)

    tspan = (0.0,70.0)
    prob = ODEProblem(twostrain_SEIR,u0,tspan,p)
    sol = solve(prob, Tsit5(), abstol = 1e-8, reltol = 1e-8, isoutofdomain=(u,p,t) -> any(x -> x < 0, u))

    return(sol)

end


function inter_annual_evolution(S_end, I1_end, I2_end, νi, νr, S_old,Z1old, Z2old, Ci, Cr, bi, br, si, sr, r, ϕ1, ϕ2, γ, sigma,stoch)

    frac_i1 = (I1_end)/(S_old)
    frac_i2 = (I2_end)/(S_old)
    Nt1 = stoch*S_end*(r + r*(si*νi + sr*νr))

    n1t1 = (1/(1 + si*νi + sr*νr))*(νi + si*(νi^2)*((bi^2)*(Ci^2) + 1) + sr*νr*νi)
    n2t1 = (1/(1 + si*νi + sr*νr))*(νr + sr*(νr^2)*((br^2)*(Cr^2) + 1) + si*νr*νi)

    #n1t1 = (1/(1 + si*νi + sr*νr))*(νi + si*(νi^2)*((bi^2)*(Ci^2) + 1) + sr*νr*νi)
    #n2t1 = (1/(1 + si*νi + sr*νr))*(νr + sr*(νr^2)*((br^2)*(Cr^2) + 1) + si*νr*νi)

    Z1t1 = ϕ1*I1_end + γ*Z1old
    Z2t1 = ϕ2*I2_end + γ*Z2old
    return frac_i1, frac_i2, Nt1, Z1t1, Z2t1, n1t1, n2t1
end

function get_tree_sp(prob_doug, n)
    wsample(["DO","GR"],[prob_doug, 1 - prob_doug],n)
end

function countmemb(itr)
    d = Dict{eltype(itr), Int}()
    for val in itr
        if isa(val, Number) && isnan(val)
            continue
        end
        d[val] = get(d, val, 0) + 1
    end
    return d
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

function dbinom(n,k,p)
    b = Binomial(k,p)
    return pdf(b,n)
end

# setting up the initial conditions
eps_val = 0.1
ν_SNPV_DO = 0.1
ν_SNPV_GR = 0.1
ν_MNPV_DO = 0.1
ν_MNPV_GR = 0.1

# C_SNPV_DO = 3.7
# C_SNPV_GR = 3.7
# C_MNPV_DO = 2.5
# C_MNPV_GR = 3.3

# C_SNPV_DO = 3.38 # model 3 no tree species
# C_SNPV_GR = 3.38 # model 3 no tree species
# C_MNPV_DO = 2.65 # model 3 no tree species
# C_MNPV_GR = 2.65 # model 3 no tree species

number_of_pops = length(Tables.matrix(unique(select(dist_df, "num"))))
#number_of_pops = 1

pop_nums = 1:number_of_pops

dist_mat = Array{Float64}(undef, number_of_pops, number_of_pops)
for i in pop_nums
    temp = select(filter(:ref=>==(i),dist_df), [:distance])
    dist_mat[:,i] = transpose(convert(Vector{Float64}, temp.distance))
end

function run_simulation(d, pdoug, eps_val, gen, phi1, phi2, νSDO, νSGR, νMDO, νMGR, CSDO, CSGR, CMDO, CMGR,rep, tree_vals, si,sr,alpha,sigma,bi,br,gamma,r)
    #tree_vals = get_tree_sp(pdoug,number_of_pops)
    #tree_vals = ["DO","GR"]
    #tree_id = DataFrame(num = pop_nums, tree_sp = tree_vals)


    if number_of_pops >2
        prop_disp = get_prop(eps_val,dist_mat)
    elseif number_of_pops == 2
        prop_disp = [0 1; 1 0]
    else
        prop_disp = Array{Float64}(undef, 1, number_of_pops)
        prop_disp[1,1] = 0
    end

    S_pop = zeros(number_of_pops, gen+1)
    Z1_pop = zeros(number_of_pops, gen+1)
    Z2_pop = zeros(number_of_pops, gen+1)
    frac_pop1 = zeros(number_of_pops, gen+1)
    frac_pop2 = zeros(number_of_pops, gen+1)
    nu1_track = zeros(number_of_pops, gen+1)
    nu2_track = zeros(number_of_pops, gen+1)

    S_pop_AD = zeros(number_of_pops, gen+1) # s population after dispersal
    Z1_pop_AD = zeros(number_of_pops, gen+1)
    Z2_pop_AD = zeros(number_of_pops, gen+1)
    nu1_AD = zeros(number_of_pops, gen+1)
    nu2_AD = zeros(number_of_pops, gen+1)

    initS = rand(Uniform(0.4,1.8),1)
    initMNPV = rand(LogNormal(-1,0.5),1)
    initSNPV = rand(LogNormal(-1,0.5),1)

    for p in 1:number_of_pops


        S_pop[p,1] = initS[1]
        Z1_pop[p,1] = initSNPV[1]
        Z2_pop[p,1] = initMNPV[1]

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

        S_pop_AD[:,t] .= transpose(S_move)
        Z1_pop_AD[:,t] .= transpose(Z1_move)
        Z2_pop_AD[:,t] .= transpose(Z2_move)
        nu1_AD[:,t] .= transpose(nu1_move)
        nu2_AD[:,t] .= transpose(nu2_move)

        en = rand(Normal(0,sigma),1)
        stoch = exp(en[1])

        for p in 1:number_of_pops

            if tree_vals[p] == "DO"
                C_SNPV = CSDO
                C_MNPV = CMDO
            else
                C_SNPV = CSGR
                C_MNPV = CMGR
            end

            S_init = S_pop_AD[p,t] #+ 1e-6
            Z1_init = Z1_pop_AD[p,t] #+ 1e-8
            Z2_init = Z2_pop_AD[p,t] #+ 1e-8
            nu1_init = nu1_AD[p,t]
            nu2_init = nu2_AD[p,t]

            output = run_ode(S_init,Z1_init,Z2_init,nu1_init,nu2_init,C_SNPV,C_MNPV,alpha)

            end_t = length(output.t)

            num_equations = size(output)[2]

            S_end = output[1,end_t]
            nu1_end = output[2,end_t]
            nu2_end = output[3,end_t]
            I1_end = output[num_equations - 1,end_t]
            I2_end = output[num_equations,end_t]

            #S_end, I1_end, I2_end, νi, νr, Z1old, Z2old, Ci, Cr, bi, br, si, sr, r, ϕ1, ϕ2, γ, sigma,stoch
            over_winter = inter_annual_evolution(S_end, I1_end, I2_end, nu1_end, nu2_end, S_init, Z1_init, Z2_init, C_SNPV, C_MNPV, bi, br, si, sr, r, phi1, phi2, gamma, sigma, stoch)
            #over_winter = inter_annual(S_end, I1_end, I2_end, lambda, phi1, phi2)

            frac_pop1[p,t] = over_winter[1]
            frac_pop2[p,t] = over_winter[2]
            S_pop[p,t] = over_winter[3]
            Z1_pop[p,t] = over_winter[4]
            Z2_pop[p,t] = over_winter[5]
            nu1_track[p,t] = over_winter[6]
            nu2_track[p,t] = over_winter[7]

        end
    end

    S_dat = DataFrame(Tables.table(transpose(S_pop)))
    rows = size(S_dat)[1]
    ncols = size(S_dat)[2]

    insertcols!(S_dat, 1, :Year => 1:rows)
    insertcols!(S_dat, 1, :Species => "S")

    Z1_dat = DataFrame(Tables.table(transpose(Z1_pop)))
    insertcols!(Z1_dat, 1, :Year => 1:rows)
    insertcols!(Z1_dat, 1, :Species => "Z1")

    Z2_dat = DataFrame(Tables.table(transpose(Z2_pop)))
    insertcols!(Z2_dat, 1, :Year => 1:rows)
    insertcols!(Z2_dat, 1, :Species => "Z2")

    fracI1_dat = DataFrame(Tables.table(transpose(frac_pop1)))
    insertcols!(fracI1_dat, 1, :Year => 1:rows)
    insertcols!(fracI1_dat, 1, :Species => "FracI1")

    fracI2_dat = DataFrame(Tables.table(transpose(frac_pop2)))
    insertcols!(fracI2_dat, 1, :Year => 1:rows)
    insertcols!(fracI2_dat, 1, :Species => "FracI2")

    nu1_dat = DataFrame(Tables.table(transpose(nu1_track)))
    insertcols!(nu1_dat, 1, :Year => 1:rows)
    insertcols!(nu1_dat, 1, :Species => "nu1")

    nu2_dat = DataFrame(Tables.table(transpose(nu2_track)))
    insertcols!(nu2_dat, 1, :Year => 1:rows)
    insertcols!(nu2_dat, 1, :Species => "nu2")

    all_dat = vcat(S_dat,Z1_dat, Z2_dat,fracI1_dat,fracI2_dat,nu1_dat,nu2_dat)
    all_long = stack(all_dat, 3:(ncols + 2))

    tree_dat = DataFrame(Tables.table(tree_vals))
    rename!(tree_dat, :Column1 => :tree_sp)
    insertcols!(tree_dat, 1, :population => 1:size(tree_dat)[1])

    tree_dat.variable = string.("Column", tree_dat.population)

    all_long = outerjoin(all_long, tree_dat, on = :variable)

    all_long.rep .= rep

    do_index = findall(x -> x =="DO", tree_vals)
    gr_index = findall(x -> x =="GR", tree_vals)

    do_snpv = frac_pop1[do_index,:]
    gr_snpv = frac_pop1[gr_index,:]

    do_mnpv = frac_pop2[do_index,:]
    gr_mnpv = frac_pop2[gr_index,:]

    total_frac = frac_pop1 + frac_pop2
    total_frac_means = mean(total_frac,dims = 1)

    f10_index = findall(x -> x >= 0.1, total_frac_means[1,:])
    f30_index = findall(x -> x >= 0.3, total_frac_means[1,:])
    fl10_index = findall(x -> x < 0.1, total_frac_means[1,:])

    s10 = mean(frac_pop1[:,f10_index])
    m10 = mean(frac_pop2[:,f10_index])

    s30 = mean(frac_pop1[:,f30_index])
    m30 = mean(frac_pop2[:,f30_index])

    sl10 = mean(frac_pop1[:,fl10_index])
    ml10 = mean(frac_pop2[:,fl10_index])

    ds10 = mean(do_snpv[:,f10_index])
    dm10 = mean(do_mnpv[:,f10_index])
    gs10 = mean(gr_snpv[:,f10_index])
    gm10 = mean(gr_mnpv[:,f10_index])

    ds30 = mean(do_snpv[:,f30_index])
    dm30 = mean(do_mnpv[:,f30_index])
    gs30 = mean(gr_snpv[:,f30_index])
    gm30 = mean(gr_mnpv[:,f30_index])

    dsl10 = mean(do_snpv[:,fl10_index])
    dml10 = mean(do_mnpv[:,fl10_index])
    gsl10 = mean(gr_snpv[:,fl10_index])
    gml10 = mean(gr_mnpv[:,fl10_index])

    means_list = [ds10,dm10,gs10,gm10,ds30,dm30,gs30,gm30,dsl10,dml10,gsl10,gml10,s10,m10,s30,m30,sl10,ml10]
    filter_list = [">=10%",">=10%",">=10%",">=10%",">=30%",">=30%",">=30%",">=30%","<10%","<10%","<10%","<10%",">=10%",">=10%",">=30%",">=30%","<10%","<10%"]
    tree_list = ["DF","DF","GR","GR","DF","DF","GR","GR","DF","DF","GR","GR","both","both","both","both","both","both"]
    morph_list = ["SNPV","MNPV","SNPV","MNPV","SNPV","MNPV","SNPV","MNPV","SNPV","MNPV","SNPV","MNPV","SNPV","MNPV","SNPV","MNPV","SNPV","MNPV"]

    summary_df = DataFrame(mean = means_list, tree = tree_list, morph = morph_list, filter = filter_list)
    S_df = DataFrame(sp = ["S","SNPV","MNPV"], mean = [mean(S_pop),mean(Z1_pop), mean(Z2_pop)],
    mean_last = [mean(S_pop[:,(gen-50):gen]), mean(Z1_pop[:,(gen-50):gen]), mean(Z2_pop[:,(gen-50):gen])],
    max_last = [maximum(S_pop[:,(gen-50):gen]), maximum(Z1_pop[:,(gen-50):gen]), maximum(Z2_pop[:,(gen-50):gen])],
    min_last = [minimum(S_pop[:,(gen-50):gen]), minimum(Z1_pop[:,(gen-50):gen]), minimum(Z2_pop[:,(gen-50):gen])],
     max = [maximum(S_pop),maximum(Z1_pop), maximum(Z2_pop)],
      min = [minimum(S_pop),minimum(Z1_pop), minimum(Z2_pop)],
      init = [initS[1],initSNPV[1], initMNPV[1]])
    S_df.sS .= si
    S_df.sM .= sr
    S_df.rep .= rep
    S_df.pdoug .= pdoug
    S_df.rho .= rho
    S_df.phi .= phi1

    summary_df.sS .= si
    summary_df.sM .= sr
    summary_df.rep .= rep
    summary_df.pdoug .= pdoug
    summary_df.rho .= rho

    total_frac = frac_pop1 + frac_pop2
    total_frac_means = mean(total_frac,dims = 1)

    frac1_means = mean(frac_pop1,dims = 1)
    frac2_means = mean(frac_pop2,dims = 1)

    snpv_frac = frac1_means ./ total_frac_means
    mnpv_frac = frac2_means ./ total_frac_means

    f10_index = findall(x -> x >= 0.1, total_frac_means[1,:])
    f30_index = findall(x -> x >= 0.3, total_frac_means[1,:])
    fl10_index = findall(x -> x < 0.1, total_frac_means[1,:])

    s10 = mean(snpv_frac[:,f10_index])
    m10 = mean(mnpv_frac[:,f10_index])

    s30 = mean(snpv_frac[:,f30_index])
    m30 = mean(mnpv_frac[:,f30_index])

    sall = mean(snpv_frac[2:(gen+1)])
    mall = mean(mnpv_frac[2:(gen+1)])

    mean_all = mean(total_frac_means)

    mean_last_S = mean(S_pop[:,(gen-50):gen])
    max_S = maximum(S_pop)
    min_last_S = minimum(S_pop[:,(gen-50):gen])

    mean_last_Z1 = mean(Z1_pop[:,(gen-50):gen])
    mean_last_Z2 = mean(Z2_pop[:,(gen-50):gen])


    if mean_last_S <= 10 && mean_last_S >= 1e-2 && max_S <= 1e2 && min_last_S >= 1e-10
        qual = "stable"
    else
        qual = "unstable"
    end

    if mean_last_Z1 <= 1e-10
        extinct_SNPV = 1
    else
        extinct_SNPV = 0
    end

    if mean_last_Z2 <= 1e-10
        extinct_MNPV = 1
    else
        extinct_MNPV = 0
    end

    small_df = ll_df

    probs_df = DataFrame()
    for i in 1:length(small_df.id)

        if qual == "stable"
            probs_list = [dbinom(small_df.MNPV[i],small_df.total[i],m10), dbinom(small_df.MNPV[i],small_df.total[i],m30),dbinom(small_df.MNPV[i],small_df.total[i],mall)]
        else
            probs_list = [exp(-30),exp(-30),exp(-30)]
        end

        prop_list = [m10,m30,mall]
        filter_list = [">=10%", ">=30%","all"]

        temp_df = DataFrame(LL = probs_list, MNPV_prop = prop_list, FI_filter = filter_list)
        temp_df.id .= small_df.id[i]

        probs_df = vcat(probs_df,temp_df)
    end

    probs_df.quality .= qual
    probs_df.sS .= si
    probs_df.sM .= sr
    probs_df.rep .= rep
    probs_df.pdoug .= pdoug
    probs_df.alpha .= alpha
    probs_df.mean_FI .= mean_all
    probs_df.extinct_S .= extinct_SNPV
    probs_df.extinct_M .= extinct_MNPV

    return summary_df, S_df, all_long, probs_df

end


b_test = [1]
sr_test = [4]
si_test = [13.75]
ndoug = 18
phi_test = [35]
#phi2_test = 20
r_test = [0.2]
gamma_test = [0.2]
sigma = 0
alpha_test = [0.5]
#is_test = 0.01:0.05:2

ll_sim_data = DataFrame()
summary_data = DataFrame()
info_data = DataFrame()
all_data = DataFrame()

for srt in sr_test
    for sit in si_test
        for gm in gamma_test
            for phi in phi_test
                for alph in alpha_test
                    for rp in 1:10

                        tree_vals = place_tree_sp(ndoug,number_of_pops)

                        sr = srt
                        si = sit
                        phi_S = phi
                        phi_M = phi
                        r = r_test[1]
                        bi = b_test[1]
                        br = b_test[1]
                        gamma = gm
                        rep = rp
                        alpha = alph

                        output_sim = run_simulation(0.2,ndoug, eps_val,200,phi_S,phi_M,ν_SNPV_DO,ν_SNPV_GR,ν_MNPV_DO,ν_MNPV_GR,C_SNPV_DO,C_SNPV_GR,C_MNPV_DO,C_MNPV_GR,rep,tree_vals,si,sr,alpha,sigma,bi,br,gamma,r)

                        summary_data = vcat(summary_data,output_sim[1])
                        info_data = vcat(info_data,output_sim[2])
                        all_data = vcat(all_data,output_sim[3])
                        ll_sim_data = vcat(ll_sim_data,output_sim[4])




                    end
                end
            end #phi
        end # gamma
    end #sit
end #srt


ll_sim_data
gd = groupby(ll_sim_data, [:sS, :sM,:id,:FI_filter])
ll_df_new = combine(gd, :LL => mean, :MNPV_prop => mean, :extinct_S => sum, :extinct_M => sum)

dir = "output_init_test/"

sim_info = "alpha"*string(myparam)*".csv"
#sim_info = "douglas_test.csv"
CSV.write(dir*"summary_"*sim_info, summary_data, header = true)
CSV.write(dir*"info_"*sim_info, info_data, header = true)
CSV.write(dir*"all_"*sim_info, all_data, header = true)
CSV.write(dir*"ll_"*sim_info, ll_sim_data, header = true)
