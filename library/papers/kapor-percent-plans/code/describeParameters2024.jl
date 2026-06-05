myround(x::Number; digits=3) = round(x,digits=digits)
myround(x; digits=4) = x

function print_results(varnames,dfmeans,dfsd,titles; f=stdout, supertitles=[])
    rvec = repeat("c",length(titles))
    println(f,"\\begin{tabular}{l$(rvec)}")
    for st in supertitles
        println(f," & "*join(st, " & ") * " \\\\")
    end
    println(f," & "*join(titles, " & ") * " \\\\")
    println(f,"\\hline")
    for r=1:length(varnames)
        println(f,varnames[r] * " & "*join(myround.(Vector( dfmeans[r,:])), " & ") * " \\\\")
        println(f,              " & ("*join(myround.(Vector( dfsd[r,:])), ") & (") * ") \\\\")
    end
    println(f, "\\hline \\hline")
    println(f,"\\end{tabular}")
end

################################################################
# 1. describe GPA production function vs. admissions 
################################################################
γadmit = θstar[1:J+nz]
γs = θstar[J+nz+1:J+nz+3]
γq!s = θstar[J+nz+4:J+nz+6]
γcfixed = θstar[J+nz+8 : J+nz+10]
γcvariable = θstar[J+nz+11 : J+nz+13]
γAppShocks = [θstar[J+nz+14]; θstar[3J+nz+19+2nx+nrc:3J+nz+20+2nx+nrc];]

df_γ_mean = DataFrame(gamma=[0;γadmit[J+1:end];0;1],
    gpa_uta=dfoutcomes.gpa_uta, persist_uta=dfoutcomes.persist_uta, stem_uta = dfoutcomes.stem_uta, persist_uta_probit=dfoutcomes.persist_uta_probit, stem_uta_probit=dfoutcomes.stem_uta_probit,
    gpa_tamu=dfoutcomes.gpa_tamu, persist_tamu=dfoutcomes.persist_tamu, stem_tamu = dfoutcomes.stem_tamu, persist_tamu_probit=dfoutcomes.persist_tamu_probit, stem_tamu_probit=dfoutcomes.stem_tamu_probit,    )


γboot = hcat( [[0; df_bootstrap_est[J+1:J+nz,Symbol("draw_$t")];0;1 ] for t=1:nbootstrap]...)
df_γ_sd = DataFrame(gamma=vec( std(γboot,dims=2,corrected=false)),)
for k in [:gpa_uta, :persist_uta, :stem_uta, :persist_uta_probit, :stem_uta_probit, 
    :gpa_tamu, :persist_tamu, :stem_tamu, :persist_tamu_probit, :stem_tamu_probit,]
    outc = hcat([df[!,k] for df in dfoutcomes_boot]...)
    df_γ_sd[!,k] = vec(std(outc,dims=2,corrected=false))
end
df_γ_p95 = DataFrame(gamma=vec( mapslices(x->quantile(x,.95),γboot,dims=2)),)
for k in [:gpa_uta, :persist_uta, :stem_uta, :persist_uta_probit, :stem_uta_probit, 
    :gpa_tamu, :persist_tamu, :stem_tamu, :persist_tamu_probit, :stem_tamu_probit,]
    outc = hcat([df[!,k] for df in dfoutcomes_boot]...)
    df_γ_p95[!,k] = vec( mapslices(x->quantile(x,.95),outc,dims=2))
end
df_γ_p5 = DataFrame(gamma=vec( mapslices(x->quantile(x,.05),γboot,dims=2)),)
for k in [:gpa_uta, :persist_uta, :stem_uta, :persist_uta_probit, :stem_uta_probit, 
    :gpa_tamu, :persist_tamu, :stem_tamu, :persist_tamu_probit, :stem_tamu_probit,]
    outc = hcat([df[!,k] for df in dfoutcomes_boot]...)
    df_γ_p5[!,k] = vec( mapslices(x->quantile(x,.05),outc,dims=2))
end
#make table of admissions and outcome parameters!
vnames_gamma = ["Constant","SAT","Class Rank","SAT Ratio","Poverty","URM","Scholarship","Caliber (q)"]
titles_gamma = ["\$ \\gamma^\\text{admit} \$","GPA", "Persist", "STEM", "Persist", "STEM", "GPA", "Persist", "STEM", "Persist", "STEM"]
 # ["\$ \\gamma^\\text{admit} \$","GPA (UTA)","Persist (UTA) LPM","STEM (UTA) LPM", "Persist (UTA) Probit", "STEM (UTA) Probit",
 #       "GPA (TAMU)","Persist (TAMU) LPM","STEM (TAMU) LPM", "Persist (TAMU) Probit", "STEM (TAMU) Probit",]
st_gamma1 =[""; repeat(["","LPM","LPM","Probit","Probit"],outer=2);]
st_gamma2 =[""; repeat(["UTA","TAMU"],inner=5);]
supertitles_gamma = [st_gamma1, st_gamma2]


inds_noprobit = [1:4; 7:9]
inds_probit = [5:6; 10:11]
open(tablespath*"/tab_gamma.tex","w") do f
    print_results(vnames_gamma,df_γ_mean[:,inds_noprobit],
    df_γ_sd[:,inds_noprobit],titles_gamma[inds_noprobit], f=f, supertitles=[stg[inds_noprobit] for stg in supertitles_gamma] )
end

open(tablespath*"/tab_gamma_probit.tex","w") do f
    print_results(vnames_gamma,df_γ_mean[:,inds_probit],
    df_γ_sd[:,inds_probit],titles_gamma[inds_probit], f=f, supertitles=[stg[inds_probit] for stg in supertitles_gamma] )
end

open(tablespath*"/tab_gamma_full.tex","w") do f
    print_results(vnames_gamma,df_γ_mean,df_γ_sd,titles_gamma, f=f, supertitles=supertitles_gamma )
end


################################################################
# 2. table of info and appcost parameters!
################################################################
vnames_info = ["Constant","Poverty","URM"]
titles_info = ["Student Info. \$ \\gamma_{s} \$", "Uncertainty \$ \\gamma_{q \\vert s} \$", "Fixed Cost \$ \\gamma^c_0 \$", "Var. Cost \$ \\gamma^c_1 \$",
    "Var. Shocks \$ \\gamma^{\\varepsilon_\\text{app}} \$"]
df_info_mean = DataFrame(
    gamma_s = θstar[J+nz+1:J+nz+3], 
    gamma_q!s = θstar[J+nz+4:J+nz+6],
    gammacfixed = θstar[J+nz+8:J+nz+10],
    gammacvariable = θstar[J+nz+11:J+nz+13],
    gammaAppShock = θstar[[J+nz+14,3J+nz+19+2nx+nrc,3J+nz+20+2nx+nrc]],
    )
df_info_sd = DataFrame(
    gamma_s = sd_θstar[J+nz+1:J+nz+3], 
    gamma_q!s = sd_θstar[J+nz+4:J+nz+6],
    gammacfixed = sd_θstar[J+nz+8:J+nz+10],
    gammacvariable = sd_θstar[J+nz+11:J+nz+13],
    gammaAppShock = sd_θstar[[J+nz+14,3J+nz+19+2nx+nrc,3J+nz+20+2nx+nrc]],
)    
open(tablespath*"/tab_info.tex","w") do f
    print_results(vnames_info,df_info_mean,df_info_sd,titles_info, f=f )
end



################################################################
# 3. information and appcost parameter figure
################################################################
fig_params = Figure(size=(1000,700))

# plot signal variances
g1 = fig_params[1,1]
pos = Axis(g1[1,1], xlabel = "HS Poverty", title = "A. Information")
xgrid = 0:0.01:1
y1 = [sqrt(log1pexp.(dot(γs,[1,x,0]))) for x in xgrid]
y1a = [sqrt(log1pexp.(dot(γs,[1,x,1]))) for x in xgrid]
y2 = [sqrt(log1pexp.(dot(γq!s,[1,x,0]))) for x in xgrid]
y2a = [sqrt(log1pexp.(dot(γq!s,[1,x,1]))) for x in xgrid]
series!(pos,xgrid,[y1 ;; y1a ;; y2 ;; y2a]',labels=["sd(s), URM=0","sd(s), URM=1","sd(q|s), URM=0","sd(q|s), URM=1"])
ylims!(low=0.0, high=1.6)
xlims!(low=-0.05,high=1.05)
#axislegend(position=:lt, nbanks=2)
g1[2,1] = Legend(g1,pos,orientation=:horizontal, nbanks=2)


# plot cost parameters
g2 = fig_params[1,2]
pos = Axis(g2[1,1],xlabel = "HS Poverty", title = "B. Application Costs")
y3 = [log1pexp.(dot(γcfixed,[1,x,0])) for x in xgrid]
y4 = [log1pexp.(dot(γcvariable,[1,x,0])) for x in xgrid]
y5 = [ std(Gumbel(0,log1pexp.(dot(γAppShocks,[1,x,0])))) for x in xgrid]
y3a = [log1pexp.(dot(γcfixed,[1,x,1])) for x in xgrid]
y4a = [log1pexp.(dot(γcvariable,[1,x,1])) for x in xgrid]
y5a = [ std(Gumbel(0,log1pexp.(dot(γAppShocks,[1,x,1])))) for x in xgrid]
series!(pos,xgrid,[y3 ;; y3a ;; y4 ;; y4a ;; y5 ;; y5a]',labels=["c0, URM=0", "c0, URM=1", "c1, URM=0", "c1, URM=1", "sd(shocks), URM=0", "sd(shocks), URM=1"])
ylims!(low=0, high=1.6)
xlims!(low=-0.05,high=1.05)

g2[2,1] = Legend(g2,pos,orientation=:horizontal,nbanks=2)

#Makie.save(tablespath*"/fig_params.png",fig_params)

# #try to draw error bars!
# function getinfoindices(θ=θstar,xgrid=0:0.01:1)
#     γs = SVector{3}(θ[J+nz+1],θ[J+nz+2],θ[J+nz+3])
#     γq!s = SVector{3}(θ[J+nz+4],θ[J+nz+5],θ[J+nz+6])
#     y1 = [sqrt(log1pexp.(dot(γs,[1,x,0]))) for x in xgrid]
#     y1a = [sqrt(log1pexp.(dot(γs,[1,x,1]))) for x in xgrid]
#     y2 = [sqrt(log1pexp.(dot(γq!s,[1,x,0]))) for x in xgrid]
#     y2a = [sqrt(log1pexp.(dot(γq!s,[1,x,1]))) for x in xgrid]
#     [y1 ;; y1a ;; y2 ;; y2a]
# end
# function getinfoerrorbounds()
#     base = getinfoindices()
#     reps = cat([getinfoindices(df_bootstrap_est[!,Symbol("draw_$t")]) for t=1:nbootstrap]..., dims=3)
#     lb = [quantile(reps[i,j,:],.05) for i=1:size(base,1), j=1:size(base,2)]
#     ub = [quantile(reps[i,j,:],.95) for i=1:size(base,1), j=1:size(base,2)]
#     return lb, ub
# end
# lb,ub = getinfoerrorbounds()
# band(0:0.01:1, lb[:,4], ub[:,4])


################################################################
# 4. admissions
################################################################
#panel 1: pr(admit | 2nd-decile, median URM and non URM applicant, poverty=p) as function of poverty
#panel 2: pr_hat(admit | z=zmedian, q=q*(z,.95)) where q* is such that pr(admit|z,q*(z,.95))=.95, as function of poverty

# pos = Axis(fig_params[1,3],xlabel= "HS Poverty", title = "Admission Chance, UT-Austin")
function getpradmit_subj(zgvec,tau,ii,θ=θstar,CC=CC)
    γs = SVector{3}(θ[J+nz+1],θ[J+nz+2],θ[J+nz+3])
    γq!s = SVector{3}(θ[J+nz+4],θ[J+nz+5],θ[J+nz+6])
    var_s = log(1 + exp(dot(γs,view(CC.zi_info,ii,:))))
    var_q!s = log(1 + exp(dot(γq!s,view(CC.zi_info,ii,:))))
    var_q = var_s + var_q!s
    coef_s!q = var_s/var_q
    var_s!q = var_s * (1-coef_s!q)
    #
    qstar = quantile(Normal(),tau) - zgvec[ii]
    psubj = 0.0
    for (_ss,ws) in zip(Qnodes,Qweights)
        ss = coef_s!q*qstar + _ss*sqrt(2)*sqrt(var_s!q)
        for (_qq,wq) in zip(Qnodes,Qweights)
            zgqj = zgvec[ii] + ss + _qq*sqrt(2)*sqrt(var_q!s)
            psubj += ccdf(Normal(),-zgqj)*ws*wq
        end
    end
    psubj
end

function get_qbins(θ=θstar, fun=mean)
    zg_eqbm = let 
        γadmit = copy( θ[1:J+nz] )
        #γadmit[1:J] .-= Δcutoff[:baseline]
        hcat( [CC.zi[ii]*γadmit for ii=1:II]... )
    end
    zg_eqbm_uta = zg_eqbm[loc_UTA,:]
    myprobs_subj = [getpradmit_subj(zg_eqbm_uta,t,ii) for t = 0.05:0.1:0.95, ii=1:II]   
    xgrid = 0:0.08:1
    hspov = ziMat[:,4]      
    qbins_ttt = hcat( [ [fun(myprobs_subj[r,CC.ttt .& (xgrid[k] .< hspov .<= xgrid[k+1])]) for r=1:size(myprobs_subj,1)] for k=1:length(xgrid)-1]...)
    qbins_rest = hcat( [ [fun(myprobs_subj[r,(.!CC.ttt) .& (xgrid[k] .< hspov .<= xgrid[k+1])]) for r=1:size(myprobs_subj,1)] for k=1:length(xgrid)-1]...)
    return qbins_ttt, qbins_rest
end

function getpradmit_true(zgvec,ii,θ=θstar,CC=CC)
    γs = SVector{3}(θ[J+nz+1],θ[J+nz+2],θ[J+nz+3])
    γq!s = SVector{3}(θ[J+nz+4],θ[J+nz+5],θ[J+nz+6])
    var_s = log(1 + exp(dot(γs,view(CC.zi_info,ii,:))))
    var_q!s = log(1 + exp(dot(γq!s,view(CC.zi_info,ii,:))))
    var_q = var_s + var_q!s
    ccdf(Normal(0,sqrt(1+var_q)), -zgvec[ii])
end
function get_obj_chance(θ=θstar)
    zg_eqbm = let 
        γadmit = copy( θ[1:J+nz] )
        hcat( [CC.zi[ii]*γadmit for ii=1:II]... )
    end
    zg_eqbm_uta = zg_eqbm[loc_UTA,:]
    myprobs = [getpradmit_true(zg_eqbm_uta,ii) for ii=1:II] 
    xgrid = 0:0.08:1
    hspov = ziMat[:,4]      
    qbins_ttt = [mean(myprobs[CC.ttt .& (xgrid[k] .< hspov .<= xgrid[k+1])]) for k=1:length(xgrid)-1]
    qbins_rest = [mean(myprobs[(.!CC.ttt) .& (xgrid[k] .< hspov .<= xgrid[k+1])]) for k=1:length(xgrid)-1]
    return qbins_ttt, qbins_rest
end

#panel 1: true chances!
q_ttt_mean, q_nottt_mean = get_obj_chance()
q_ttt_boot = hcat( [ get_obj_chance(df_bootstrap_est[!,Symbol("draw_$t")])[1] for t=1:nbootstrap]...)
lowerrors1_ttt = q_ttt_mean .- [quantile(q_ttt_boot[t,:],.025) for t=1:size(q_ttt_boot,1)]
higherrors1_ttt = [quantile(q_ttt_boot[t,:],.975) for t=1:size(q_ttt_boot,1)] .- q_ttt_mean
q_nottt_boot = hcat( [ get_obj_chance(df_bootstrap_est[!,Symbol("draw_$t")])[2] for t=1:nbootstrap]...)
lowerrors1_nottt = q_nottt_mean .- [quantile(q_nottt_boot[t,:],.025) for t=1:size(q_nottt_boot,1)]
higherrors1_nottt = [quantile(q_nottt_boot[t,:],.975) for t=1:size(q_nottt_boot,1)] .- q_nottt_mean

#panel 2: subj chances
qb_ttt_mean, qb_nottt_mean = get_qbins()
qb_ttt_boot = hcat( [get_qbins(df_bootstrap_est[!,Symbol("draw_$t")])[1][10,:] for t=1:nbootstrap]...)
lowerrors_ttt = qb_ttt_mean[10,:] .- [quantile(qb_ttt_boot[t,:],.025) for t=1:size(qb_ttt_boot,1)]
higherrors_ttt = [quantile(qb_ttt_boot[t,:],.975) for t=1:size(qb_ttt_boot,1)] .- qb_ttt_mean[10,:]
qb_nottt_boot = hcat( [get_qbins(df_bootstrap_est[!,Symbol("draw_$t")])[2][10,:] for t=1:nbootstrap]...)
lowerrors_nottt = qb_nottt_mean[10,:] .- [quantile(qb_nottt_boot[t,:],.025) for t=1:size(qb_nottt_boot,1)]
higherrors_nottt = [quantile(qb_nottt_boot[t,:],.975) for t=1:size(qb_nottt_boot,1)] .- qb_nottt_mean[10,:]

#make the figure!
#fig_admissions = Figure(size=(1000,500))
#pos = Axis(fig_admissions[1,1], xlabel = "HS Poverty", title = "True Admission Chance, UT-Austin")
xgr = 0:0.08:1
pos = Axis(fig_params[2,1], xlabel = "HS Poverty", title = "C. True Admission Chance, UT-Austin")
scatter!(pos,xgr[1:end-1],q_ttt_mean)
errorbars!(pos,xgr[1:end-1],q_ttt_mean,lowerrors1_ttt,higherrors1_ttt)
scatter!(pos,xgr[1:end-1],q_nottt_mean)
errorbars!(pos,xgr[1:end-1],q_nottt_mean,lowerrors1_nottt,higherrors1_nottt)
ylims!(low=.35,high=1.1)
xlims!(low=-0.05,high=1.05)


#panel 2: subj chances | truth = .95
#pos = Axis(fig_admissions[1,2], xlabel = "HS Poverty", title = "Subj. Admission Chance, UT-Austin | True Chance = .95")
pos = Axis(fig_params[2,2], xlabel = "HS Poverty", title = "D. Subj. Admission Chance, UT-Austin | True Chance = .95")
scatter!(pos,xgr[1:end-1],qb_ttt_mean[10,:])
errorbars!(pos,xgr[1:end-1],qb_ttt_mean[10,:],lowerrors_ttt,higherrors_ttt)
scatter!(pos,xgr[1:end-1],qb_nottt_mean[10,:])
errorbars!(pos,xgr[1:end-1],qb_nottt_mean[10,:],lowerrors_nottt,higherrors_nottt)
ylims!(low=.35,high=1.1)
xlims!(low=-0.05,high=1.05)

colors = Makie.wong_colors()
labels = ["Top 10%","Other 90%"]
el = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

#reuse (el, labels) from fig_prefs_aware
#Legend(fig_admissions[2,1:2], el ,labels, orientation=:horizontal)
Legend(fig_params[3,1:2], el ,labels, orientation=:horizontal)
#Makie.save(tablespath*"/fig_admissions.png",fig_admissions)
Makie.save(tablespath*"/fig_params.png",fig_params)




#let qstar(z,tau) be such that pr(zgamma + q + nu > pi_uta |z,q=qstar_95(z,tau)) = tau. 
#then compute \int pr(zgamma + q + nu > pi_uta | z, q) f_{q \vert s}(q|s) f_{s\vert q}(s | qstar_95(z,tau)) ds
#this is the avg. subjective belief that i will be admitted to j, when his true chance is in fact 95%.

################################################################
# 5. awareness and preference figure
################################################################
function getpraware_simple(θ,ii,mm,CC=CC)
    grid = CC.grid[ii]
    βx_aware = view(θ,J+nz+17+nx+nrc:J+nz+16+2nx+nrc)
    βy_aware = θ[J+nz+17+2nx+nrc]
    σ_βi0_aware = exp(θ[J+nz+18+2nx+nrc])
    pr_aware = zeros(J)
    for (qq,ww) in zip(Qnodes,Qweights)
        βi0aware = qq*sqrt(2)*σ_βi0_aware
        x_terms = CC.xij[ii]*βx_aware
        pr_aware += ww .* logistic.( x_terms .+ βy_aware*grid.y[mm] .+ βi0aware)
    end
    pr_aware
end
prAwareIJ = zeros(J,II)
for ii=1:II
    for mm=1:ndraws
        prAwareIJ[:,ii] .+= getpraware_simple(θstar,ii,mm)./ndraws
    end
end

#misspecified "ignore selection" version
prAwareIJ_n = zeros(J,II)
for ii=1:II
    for mm=1:ndraws
        prAwareIJ_n[:,ii] .+= getpraware_simple(θ0,ii,mm)./ndraws
    end
end

#1st-preference shares with and without aid
function getu_simple(θ,ii,mm,aware,CC=CC)
    βp = SVector{2}(θ[J+nz+15],θ[J+nz+16])
    βx = view(θ,J+nz+17 : J+nz+16+nx) #deterministic coefficients
    σ_rc = exp.( view(θ,J+nz+17+nx: J+nz+16+nx+nrc))
    log_σ_e = view(θ,2J+nz+19+2nx+nrc : 3J+nz+18+2nx+nrc)
    coef_aidamount = view(θ,J+nz+18+2nx+nrc+1:2J+nz+18+2nx+nrc)
    grid = CC.grid[ii]
    aid_m = max.(0, (CC.listprice .- grid.efc[mm]).* exp.(coef_aidamount)) .* aware
    netprice = CC.listprice .- aid_m
    pricecoef = -log(1+exp(βp[1] + βp[2]*grid.y[mm]))
    rc = σ_rc .* view(grid.rc,:,mm)
    σ_e = exp.(log_σ_e)
    u = CC.xij[ii]*βx .+ CC.xrc*rc
    u .+= netprice .* pricecoef
    u .+= σ_e .* view(grid.ep,:,mm)
    u
end
function enrollmentprob_simple(ui,θ,portf_B,CC=CC)
    λ = 1/(1+exp(-θ[J+nz+7]))
    umax = maximum(ui)
    expλu = exp.((ui .- umax)./λ)
    inner_terms =portf_B .* expλu
    sum_inner_terms = sum(inner_terms)
    _G = sum_inner_terms^λ + exp(-umax)
    inner_terms .* ( sum_inner_terms^(λ - 1) / _G )
end

prIJ0= zeros(J,II)
prIJ1= zeros(J,II)
prIJ1all= zeros(J,II)
for ii=1:II
    mod(ii,1000)==1 && println("getting choice probs, ii=$ii")
    choiceset = [1,1,1,1,0,1,1]
    for mm=1:ndraws
        ui = getu_simple(θstar,ii,mm,falses(7))
        prj = enrollmentprob_simple(ui,θstar,choiceset)
        prIJ0[:,ii] += prj./ndraws
        #
        ui = getu_simple(θstar,ii,mm,trues(7))
        prj = enrollmentprob_simple(ui,θstar,choiceset)
        prIJ1[:,ii] += prj./ndraws
        prIJ1all[:,ii] += enrollmentprob_simple(ui,θstar,ones(Int,J))./ndraws
    end
end

prIJ0_n= zeros(J,II)
prIJ1_n = zeros(J,II)
prIJ1all_n = zeros(J,II)
for ii=1:II
    mod(ii,1000)==1 && println("getting choice probs, ii=$ii")
    choiceset = [1,1,1,1,0,1,1]
    for mm=1:ndraws
        ui = getu_simple(θ0,ii,mm,falses(7))
        prj = enrollmentprob_simple(ui,θ0,choiceset)
        prIJ0_n[:,ii] += prj./ndraws
        #
        ui = getu_simple(θ0,ii,mm,trues(7))
        prj = enrollmentprob_simple(ui,θ0,choiceset)
        prIJ1_n[:,ii] += prj./ndraws
        prIJ1all_n[:,ii] += enrollmentprob_simple(ui,θ0,ones(Int,J))./ndraws
    end
end

### make some bar charts: shares; shares dropping selective private; pr(aware)
fig_prefs_aware = Figure(size=(1000,475))
category = repeat(1:J,inner=4)
grp = repeat(1:4,outer=J)

pos = Axis(fig_prefs_aware[1,1], title = "A. First Choice If Aware of Aid", xticks=(1:J,["($j)" for j=1:J]))
elements = vec( [mean(prIJ1all[:,CC.ttt],dims=2);; mean(prIJ1all[:,.!CC.ttt],dims=2) ;; mean(prIJ1all_n[:,CC.ttt],dims=2);; mean(prIJ1all_n[:,.!CC.ttt],dims=2)]')
barplot!(pos,category,elements,dodge=grp,color=colors[grp])

pos = Axis(fig_prefs_aware[1,2], title = "B. Share Aware of Aid", xticks=(1:J,["($j)" for j=1:J]))
elements = vec( [mean(prAwareIJ[:,CC.ttt],dims=2);; mean(prAwareIJ[:,.!CC.ttt],dims=2) ;; mean(prAwareIJ_n[:,CC.ttt],dims=2);; mean(prAwareIJ_n[:,.!CC.ttt],dims=2)]')
barplot!(pos,category,elements,dodge=grp,color=colors[grp])

pos = Axis(fig_prefs_aware[2,1], title = "C. First Choice, Excluding Elite Private, If Aware of Aid", xticks=(1:J,["($j)" for j=1:J]))
elements = vec( [mean(prIJ1[:,CC.ttt],dims=2);; mean(prIJ1[:,.!CC.ttt],dims=2) ;; mean(prIJ1_n[:,CC.ttt],dims=2);; mean(prIJ1_n[:,.!CC.ttt],dims=2)]')
barplot!(pos,category,elements,dodge=grp,color=colors[grp])
ylims!(low=-0.01, high=0.3)

pos = Axis(fig_prefs_aware[2,2], title = "D. First Choice, Excluding Elite Private, If Not Aware of Aid", xticks=(1:J,["($j)" for j=1:J]))
elements = vec( [mean(prIJ0[:,CC.ttt],dims=2);; mean(prIJ0[:,.!CC.ttt],dims=2) ;; mean(prIJ0_n[:,CC.ttt],dims=2);; mean(prIJ0_n[:,.!CC.ttt],dims=2)]')
barplot!(pos,category,elements,dodge=grp,color=colors[grp])
ylims!(low=-0.01, high=0.3)


labels2 = ["Top 10%","Other 90%","Top 10% (Ignore Selection)", "Other 90% (Ignore Selection)"]
el2 = [PolyElement(polycolor = colors[i]) for i in 1:length(labels2)]
Legend(fig_prefs_aware[3,1:2], el2 ,labels2, orientation=:horizontal)

Makie.save(tablespath*"/fig_prefs_aware.png",fig_prefs_aware)












### make some bar charts: shares; shares dropping selective private; pr(aware)
fig_prefs_aware_simple = Figure(size=(1000,475))
category = repeat(1:J,inner=2)
grp = repeat(1:2,outer=J)

pos = Axis(fig_prefs_aware_simple[1,1], title = "A. First Choice If Aware of Aid", xticks=(1:J,["($j)" for j=1:J]))
elements = vec( [mean(prIJ1all[:,CC.ttt],dims=2);; mean(prIJ1all[:,.!CC.ttt],dims=2)]')
barplot!(pos,category,elements,dodge=grp,color=colors[grp])

pos = Axis(fig_prefs_aware_simple[1,2], title = "B. Share Aware of Aid", xticks=(1:J,["($j)" for j=1:J]))
elements = vec( [mean(prAwareIJ[:,CC.ttt],dims=2);; mean(prAwareIJ[:,.!CC.ttt],dims=2)]')
barplot!(pos,category,elements,dodge=grp,color=colors[grp])

pos = Axis(fig_prefs_aware_simple[2,1], title = "C. First Choice, Excluding Elite Private, If Aware of Aid", xticks=(1:J,["($j)" for j=1:J]))
elements = vec( [mean(prIJ1[:,CC.ttt],dims=2);; mean(prIJ1[:,.!CC.ttt],dims=2)]')
barplot!(pos,category,elements,dodge=grp,color=colors[grp])
ylims!(low=-0.01, high=0.3)

pos = Axis(fig_prefs_aware_simple[2,2], title = "D. First Choice, Excluding Elite Private, If Not Aware of Aid", xticks=(1:J,["($j)" for j=1:J]))
elements = vec( [mean(prIJ0[:,CC.ttt],dims=2);; mean(prIJ0[:,.!CC.ttt],dims=2)]')
barplot!(pos,category,elements,dodge=grp,color=colors[grp])
ylims!(low=-0.01, high=0.3)


labels1 = ["Top 10%","Other 90%"]
el1 = [PolyElement(polycolor = colors[i]) for i in 1:length(labels1)]
Legend(fig_prefs_aware_simple[3,1:2], el1 ,labels1, orientation=:horizontal)

Makie.save(tablespath*"/fig_prefs_aware_simple.png",fig_prefs_aware_simple)















#numbers for paper 
mean(prAwareIJ[:,CC.ttt] ) #81%
mean(prAwareIJ[:,.!CC.ttt]) #62%
sum( mean(prIJ1[:,CC.ttt],dims=2)) #92.5% XXX GET SD XXX
sum( mean(prIJ1[:,.!CC.ttt],dims=2)) #92%
sum( mean(prIJ0[:,CC.ttt],dims=2)) #90.4% XXX GET SD XXX
sum( mean(prIJ0[:,.!CC.ttt],dims=2)) #89.5%
sum( mean(prIJ1,dims=2)) #92.1%
sum( mean(prIJ0,dims=2)) #89.8%
sum( mean(prIJ1,dims=2)) - sum( mean(prIJ0,dims=2)) #2.3%





#####################
# change in pr(apply), pr(admit), pr(enroll), pr(persist) by GPA 
#where "pulled from"?


################################################################
# 7. extra stuff
################################################################
# #appendix figure: total 4-year enrollment
# newcols = 1:J
# enroll_top10_all = OrderedDict(k=>100*sum(v[newcols,CC.ttt]) for (k,v) in shares)
# enroll_rest_all = OrderedDict(k=>100*sum(v[newcols,.!(CC.ttt)]) for (k,v) in shares)
# enroll_urm_all = OrderedDict(k=>100*sum(v[newcols,urm]) for (k,v) in shares)
# enroll_theil_all = OrderedDict(k=>theilIndex(v,newcols) for (k,v) in shares)
# plots_a1_all = []
# for (res,name) in zip([enroll_top10_all,enroll_rest_all,enroll_urm_all,enroll_theil_all],
#     ["Top Decile", "Not Top Decile", "URM", "Theil Index:"])
#     ys = collect(values(res))
#     plt = waterfall(1:4,[ys[1]; ys[2:end].-ys[1:end-1]],axis = (xticks = (1:4, xlabels), title="$name Total Enrollment"),)
#     push!(plots_a1_all,plt)
# end



#experimental stuff: plot costs
xgrid = 0:0.01:1
γcfixed = θstar[J+nz+8:J+nz+10]
γcvariable = θstar[J+nz+11:J+nz+13]
y1 = [log1pexp.(dot(γcfixed,[1,x,0])) for x in xgrid]
y2 = [log1pexp.(dot(γcvariable,[1,x,0])) for x in xgrid]
plt,ax,sp = series(xgrid,[y1 ;; y2]',labels=["c0 (fixed cost)","c1 (variable cost)"])
axislegend(ax)
plt

#plot mean utilities
βx = θstar[J+nz+17 : J+nz+16+nx]
meanu = sum( popwt[k] .* (CC.xij[k]*βx) for k=1:II) ./ sum(popwt)
barplot(meanu)

#plot admission cutoffs
γadmit = θstar[1:J+nz]
barplot(γadmit[1:J])

#attendance prob as function of test score!
satr = [CC.zi[ii][1,J+1] for ii=1:II]
classrank = [CC.zi[ii][1,J+2] for ii=1:II]

#E(q|enroll) binscatter
myranks = round.( 1 .- sort(unique(classrank)), digits=1)
Eq_uta = [mean(Eq[:baseline][7,findall(classrank.==k)], Weights(shares[:baseline][7,findall(classrank.==k)])) for k = sort(unique(classrank))]
Eq_tamu = [mean(Eq[:baseline][6,findall(classrank.==k)], Weights(shares[:baseline][6,findall(classrank.==k)])) for k = sort(unique(classrank))]
scatter(myranks,Eq_uta)
scatter(myranks,Eq_tamu)

Ev_uta = [mean(Ev[:baseline][7,findall(classrank.==k)], Weights(shares[:baseline][7,findall(classrank.==k)])) for k = sort(unique(classrank))]
scatter(myranks,Ev_uta)



function mybins(series,wseries,runningvar,_inds_=true)
    by = sort(unique(runningvar))
    inds = [findall(_inds_ .& (runningvar.==k)) for k in by]
    [mean(series[ix],Weights(wseries[ix])) for ix in inds]
end
Eq_uta_lowSAT = mybins(Eq[:baseline][7,:],shares[:baseline][7,:],classrank, satr .< median(satr))
Eq_uta_highSAT = mybins(Eq[:baseline][7,:],shares[:baseline][7,:],classrank, satr .>= median(satr))



### big param tables!
function report_θ(df=df_bootstrap_est)
    name_x = ["Dist < 25","Distance","URM","URM X UTA", "URM X TAMU", "Poverty","Poverty X UTA", "Poverty X TAMU", "LOS X UTA", "Century X TAMU", "SAT", "Class Rank", "SAT / HS Mean SAT", "SAT X TAMU", "SAT X UTA", "SAT X Private"]
    name_z = ["SAT","Class Rank", "SAT/HS mean SAT","Poverty","URM"]
    name_J = ["In-State Public","Private","Out-of-State Public","Relig.","Elite Private","Texas A\\&M", "UT Austin"]
    name_rc = ["Distance","S/F Ratio","UTA vs TAMU"]
    name_zinfo = ["Const","Poverty","URM"]
    #
    paramnames = [
        ["\$\\underline\\pi\$ Eqbm. ($(z))" for z in name_J];
        ["\$\\gamma\$ ($z)" for z in name_z];
        ["\$\\gamma^s\$ ($z)" for z in name_zinfo];
        ["\$\\gamma^{q|s}\$ ($z)" for z in name_zinfo];
        "\$\\lambda\$ (Matric. shock scale)";
        ["\$\\gamma^{fixed}\$ ($z)" for z in name_zinfo];
        ["\$\\gamma^{var}\$ ($z)" for z in name_zinfo];
        ["\$\\gamma^{shock}\$ ($z)" for z in name_zinfo[1:1]];
        ["\$\\beta^{p}\$ ($z)" for z in ["Const","income"]];
        ["\$\\beta^{(w,x,z)}\$ ($x)" for x in name_J];
        ["\$\\beta^{(w,x,z)}\$ ($x)" for x in name_x];
        ["\$\\log(\\sigma^{rc})\$ ($x)" for x in name_rc];
        ["\$\\beta^{aware}\$ ($x)" for x in name_J];
        ["\$\\beta^{aware}\$ ($x)" for x in name_x];
        "\$\\beta^{aware}\$ (income)";
        "\$\\log(\\sigma) (\\beta_i^0)\$";
        ["\$\\alpha^{aid}\$ ($z)" for z in name_J];
        ["\$\\log(\\sigma^{e})\$ ($x)" for x in name_J];
        ["\$\\gamma^{shock}\$ ($z)" for z in name_zinfo[2:3]];
    ]
    #
    θ = df.thetastar
    θboot = Array(df[:,2:end])
    se = vec(std(θboot,dims=2))
    p95 = [quantile(view(θboot,m,:),.95) for m=1:length(θ)]
    p5 = [quantile(view(θboot,m,:),.05) for m=1:length(θ)]
    df_out = DataFrame(ind=1:length(θ),paramnames=paramnames,estimate=θ,se=se,p5=p5,p95=p95)
    param_inds = [ (1:J+nz+14); (3J+nz+19+2nx+nrc:3J+nz+20+2nx+nrc);  (J+nz+15:3J+nz+18+2nx+nrc); ]
    df_out = df_out[param_inds,:]
    df_out.ind .= 1:length(θ)
    df_out
end

myround(x::Int; digits=2) =x
function print_params(df,f=stdout)
    rvec = repeat("r",4)
    println(f,"\\begin{tabular}{ll$(rvec)}")
    println(f," & Parameter & Estimate & SE & p5 & p95 \\\\")
    for row in eachrow(df)
        println(f, join(myround.(Vector(row)), " & ") * " \\\\")
    end
    println(f, "\\hline \\hline")
    println(f,"\\end{tabular}")
end
    
dfparams = report_θ()
dfparams_1 = dfparams[1:32,:]
dfparams_2 = dfparams[33:64,:]
dfparams_3 = dfparams[65:end,:]

for (t,df) in enumerate([dfparams_1,dfparams_2,dfparams_3])  
    open(tablespath*"/tab_params_$t.tex","w") do f
        print_params(df,f)
    end
end


###cutoffs 
Δcutoff_boot = []
for t=1:nbootstrap
    push!(Δcutoff_boot, JLD2.load("$filepath/bootstrap_$(t)_step3.jld2","Δcutoff"))
end

function get_df_cutoff(c)
    mynames = ["TTP", "Aid Awareness Only", "Aid Awareness + TTP", "Top 5\\% zg", "Top 10\\% zg", "Top 15\\% zg", "Top 20\\% zg"]
    myfun = ind -> ([0; c[:allaware][ind]; c[:TTT_allaware][ind]; c[:by_zg_95][ind]; c[:by_zg_90][ind]; c[:by_zg_85][ind]; c[:by_zg_80][ind];] .- c[:baseline][ind])
    DataFrame(Policy = mynames, UTA = myfun(7), TAMU = myfun(6), otherpublic = myfun(1))
end

dfcutoff_point = get_df_cutoff(Δcutoff)
dfcutoff_sd = let 
    dfcs = [get_df_cutoff(c) for c in Δcutoff_boot]
    Dict(k=> std([df[!,k] for df in dfcs]) for k in [:UTA,:TAMU,:otherpublic]) |> DataFrame
end

function print_cutoffs(dfcutoff_point=dfcutoff_point,dfcutoff_sd=dfcutoff_sd; f=stdout)
    println(f,"\\begin{tabular}{lrrr}")
    println(f,"Policy & UT Austin & Texas A\\&M & Other Public \\\\")
    for r = 1:nrow(dfcutoff_point)
        println(f, join(myround.(Vector(dfcutoff_point[r,:])), " & ") * " \\\\")
        println(f, " & ("*join(myround.(Vector(dfcutoff_sd[r,:]),digits=2), ") & (") * ") \\\\ \\noalign{\\vskip 2mm} ")
    end
    println(f, "\\hline \\hline")
    println(f,"\\end{tabular}")
end

open(tablespath*"/tab_cutoffs.tex","w") do f
    print_cutoffs(f=f)
end