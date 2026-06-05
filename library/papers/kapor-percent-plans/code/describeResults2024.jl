


########
# describe counterfactuals
########

##################################
#= construct sum stats:
    #Enrollment (top 10%);
    #Enrollment (rest)
    #Enriollment (% URM)
    #theil index of high schools
    #GPA 
    #Persist
    #STEM major
    #Caliber
    #program payoffs
=# ###############################
mycols = 6:7
include("describeResultsFunctions2024.jl")

sumstats = summarize_cfs()

#try w/ bootstrap draws
sumstats_boot = []
dfr_boot = []
popwts_boot = []
inds_boot = []
for t=1:nbootstrap
    println("summarizing bootstrap rep $t")
    stuff = JLD2.load("$filepath/bootstrap_$t.jld2")
    stuff_step3 = JLD2.load("$filepath/bootstrap_$(t)_step3.jld2")
    for k in keys(stuff_step3)
        @assert k in keys(stuff)
        stuff[k] = stuff_step3[k]
    end
    push!(inds_boot,stuff["inds"])
    θb = df_bootstrap_est[:,Symbol("draw_$t")]
    newpopwt = (popwt .* df_bootstrap_wts[:,Symbol("draw_$t")])[stuff["inds"]]
    newpopwt ./= sum(newpopwt)
    push!(popwts_boot,newpopwt)
    sumstats_t = summarize_cfs(stuff["CC"],stuff["inds"],θb,stuff["θoutcome"],  
        stuff["Eq"],stuff["Ev"],stuff["Eshock"],stuff["Δcutoff"],stuff["shares"], newpopwt,)
    push!(sumstats_boot, sumstats_t)
    zg_uta_t = let
        γadmit = θb[1:J+nz]
        γadmit[1:J] .-= stuff["Δcutoff"][:baseline]
        zg = hcat( [stuff["CC"].zi[ii]*γadmit for ii=1:length(stuff["inds"])]... )
        zg[loc_UTA,:]
    end
    push!(dfr_boot, make_robustness_table(sumstats_t,stuff["CC"],zg_uta_t,newpopwt))
end

sd_sumstats = OrderedDict(k=> OrderedDict(k2=> std( [st[k][k2] for st in sumstats_boot] ) for k2 in keys(sumstats_boot[1][k])) for k in keys(sumstats_boot[1]))

xlabels = ["baseline","direct","information","equilibrium"]
ys = collect( values(sumstats[:enroll_theil]))


function compute_x(x, dodge, width=1, gap=0.2, dodge_gap=0.03) #for adding error bars to grouped barcharts; makie default options
    scale_width(dodge_gap, n_dodge) = (1 - (n_dodge - 1) * dodge_gap) / n_dodge
    function shift_dodge(i, dodge_width, dodge_gap)
        (dodge_width - 1) / 2 + (i - 1) * (dodge_width + dodge_gap)
    end
    width *= 1 - gap
    n_dodge = maximum(dodge)
    dodge_width = scale_width(dodge_gap, n_dodge)
    shifts = shift_dodge.(dodge, dodge_width, dodge_gap)
    return x .+ width .* shifts
end

#cf fig 1: flagship enrollment, combine "base" and "with info
urm = ziRaw[:,:minority] .==1
highpoverty = ziRaw[:,:econdisadv] .>= quantile(ziRaw.econdisadv,.75)
affluent = ziRaw[:,:econdisadv] .<= quantile(ziRaw.econdisadv,.25)
denoms_fig1 = [sum(popwt[CC.ttt]); sum(popwt[urm]); sum(popwt[highpoverty]); sum(popwt[affluent]); 1.0] 
stats_fig1 = [:enroll_top10,:enroll_urm,:enroll_highpoverty,:enroll_affluent,:enroll_theil]
basevals = [round(sumstats[_res][:baseline]/denoms_fig1[m],digits=2) for (m,_res) in enumerate(stats_fig1)]
names_fig1 = ["A. Flagship Enrollment: Top Decile (% Change, Baseline=$(basevals[1]))", "B. Flagship Enrollment: URM (% Change, Baseline=$(basevals[2]))", "C. Flagship Enrollment: High-Poverty HS (Baseline=$(basevals[3]))", "D. Flagship Enrollment: Affluent HS (Baseline=$(basevals[4]))", "E. Flagship Enrollment: HS Concentration (Theil Index, Baseline=$(basevals[5]))"]
xlab2 = ["Baseline","Aid Info: Direct","Aid Info: Eqbm.","Mechanical","Information","Eqbm."]
colors = Makie.wong_colors()
fig_1 = Figure(size=(1000,1000))
for (nn,_res,name) in zip(1:5, stats_fig1, names_fig1)
    res = sumstats[_res]
    ys = [res[:baseline], res[:direct], res[:info], res[:eqbm]]
    yvals_1 = [100*(ys[k]-ys[k-1])./ys[1] for k=2:length(ys)]
    yvals_1_boot = zeros(length(yvals_1),nbootstrap)
    for t=1:nbootstrap
        res_b = sumstats_boot[t][_res]
        ys = [res_b[:baseline], res_b[:direct], res_b[:info], res_b[:eqbm]]
        yvals_1_b = [100*(ys[k]-ys[k-1])./ys[1] for k=2:length(ys)]
        yvals_1_boot[:,t] .= yvals_1_b
    end
    yv1_low = yvals_1 .- [quantile(yvals_1_boot[r,:],.025) for r=1:size(yvals_1_boot,1)]
    yv1_hi = [quantile(yvals_1_boot[r,:],.975) for r=1:size(yvals_1_boot,1)] .- yvals_1
    #
    ys_2 = [res[:baseline], res[:nottt_allaware], res[:base_allaware], res[:direct_allaware], res[:info_allaware], res[:eqbm_allaware]]
    yvals_2 = [100*(ys_2[k]-ys_2[k-1])./ys_2[1] for k=2:length(ys_2)]
    #
    yvals_2_boot = zeros(length(yvals_2),nbootstrap)
    for t=1:nbootstrap
        res_b = sumstats_boot[t][_res]
        ys_2 = [res_b[:baseline], res_b[:nottt_allaware], res_b[:base_allaware], res_b[:direct_allaware], res_b[:info_allaware], res_b[:eqbm_allaware]]
        yvals_2_b = [100*(ys_2[k]-ys_2[k-1])./ys_2[1] for k=2:length(ys_2)]
        yvals_2_boot[:,t] .= yvals_2_b
    end
    yv2_low =  yvals_2 .- [quantile(yvals_2_boot[r,:],.025) for r=1:size(yvals_2_boot,1)] 
    yv2_hi = [quantile(yvals_2_boot[r,:],.975) for r=1:size(yvals_2_boot,1)] .- yvals_2
    #
    mykeys = [:baseline_by_zg_80,:direct_by_zg_80,:info_by_zg_80,:by_zg_80]
    ys = [res[k] for k in mykeys]
    yvals_3 = [100*(ys[k]-ys[k-1])./ys[1] for k=2:length(ys)]
    yvals_3_boot = zeros(length(yvals_3),nbootstrap)
    for t=1:nbootstrap
        res_b = sumstats_boot[t][_res]
        ys = [res_b[k] for k in mykeys]
        yvals_3_b =  [100*(ys[k]-ys[k-1])./ys[1] for k=2:length(ys)]
        yvals_3_boot[:,t] .= yvals_3_b
    end
    yv3_low = yvals_3 .- [quantile(yvals_3_boot[r,:],.025) for r=1:size(yvals_3_boot,1)]
    yv3_hi = [quantile(yvals_3_boot[r,:],.975) for r=1:size(yvals_3_boot,1)] .- yvals_3
    #
    mykeys = [:baseline_by_zg_80,:direct_by_zg_80,:info_by_zg_80,:by_zg_80]
    ys = [res[k] for k in mykeys]
    yvals_4 = [100*(ys[k]-ys[k-1])./ys[1] for k=2:length(ys)]
    yvals_4_boot = zeros(length(yvals_3),nbootstrap)
    for t=1:nbootstrap
        res_b = sumstats_boot[t][_res]
        ys = [res_b[k] for k in mykeys]
        yvals_3_b =  [100*(ys[k]-ys[k-1])./ys[1] for k=2:length(ys)]
        yvals_3_boot[:,t] .= yvals_3_b
    end
    yv3_low = yvals_3 .- [quantile(yvals_3_boot[r,:],.025) for r=1:size(yvals_3_boot,1)]
    yv3_hi = [quantile(yvals_3_boot[r,:],.975) for r=1:size(yvals_3_boot,1)] .- yvals_3
    #
    yvals = [yvals_1; yvals_2; yvals_3]
    yv_low = [yv1_low; yv2_low; yv3_low]
    yv_hi = [yv1_hi; yv2_hi; yv3_hi]
    ll = min( minimum( cumsum(yvals_1)), minimum(cumsum(yvals_2)), -15) * 1.5
    ul = max( maximum( cumsum(yvals_1)), maximum(cumsum(yvals_2)), 10) * 1.5
    nn==2 && (ul=40)
    pos = fig_1[nn,1]
    ngroup = 3
    group = [fill(1,3); fill(2,5); fill(3,3)] #repeat(1:ngroup, inner=length(yvals_2))
    x = [3,4,5,1,2,3,4,5,3,4,5] #repeat(1:length(yvals_2), outer=ngroup)
    #xlab = repeat(xlab2[2:end],outer=ngroup)
    xlab = xlab2[2:end][x]
    #
    _ax,_wf = waterfall(pos, x, yvals ,axis = (xticks = (x, xlab), title="$(name)", limits = (nothing,nothing,ll,ul)),
    show_direction=false,
    # color_over_bar=:black,
    # color_over_background=:black,
    # bar_labels=[(y>0 ? "+"*string( round(y,digits=2))*"%" : 
    #              y<0 ? string( round(y,digits=2))*"%" :
    #              "") for y in yvals],
    dodge=group, color=colors[group], stack=x )
    bar_labels=[(y>0 ? "+"*string( round(y,digits=2))*"%" : 
                 y<0 ? string( round(y,digits=2))*"%" :
                 "") for y in yvals]
    xerr = compute_x(x,group)
    errorbars!(xerr, [cumsum(yvals_1); cumsum(yvals_2); cumsum(yvals_3)], yv_low, yv_hi,
    whiskerwidth=5, color=:gray )
    #hack to get waterfall plot to display bar labels
    wf_args = _wf.plots[1].args
    wf_attr = _wf.plots[1].attributes
    delete!(_ax,_wf)
    barplot!(_ax,wf_args...; wf_attr..., bar_labels=bar_labels,
    dodge=group, color=colors[group], stack=x,
    color_over_background=:black, color_over_bar=:black,)
end
Legend(fig_1[6,1], [PolyElement(polycolor=colors[i]) for i=1:3], ["TTP","TTT+Aid Awareness","Alternative Weights, Top 20%"], orientation=:horizontal)

Makie.save(tablespath*"/fig_cf_1.png",fig_1)


# #cf fig 1 part b: academic outcomes (barchart "pulled in / pulled out" version)
xlab3 = ["Baseline"; "Aid Info"; "Mechanical"; "Information"; "Eqbm"]
fig_2 = Figure(size=(1000,750))
for (nn,_res,name) in zip(1:4,[:gpa_flagship,:persist_flagship,:stem_flagship,:E_pi_flagship,],
    ["GPA", "1(Persist)", "1(STEM)", "Program Payoffs",])
    rows = 6:7
    res = sumstats[_res]
    ks = [:baseline,:direct,:info,:eqbm]
    yvals_1 = [res[ks[1]]; [(res[ks[k]]*sum(shares[ks[k]][rows,:])-res[ks[k-1]]*sum(shares[ks[k-1]][rows,:]))/sum(shares[ks[k]][rows,:]-shares[ks[k-1]][rows,:]) for k=2:length(ks)];]
    #
    ks_2 = [:baseline,:nottt_allaware,:base_allaware,:direct_allaware,:info_allaware,:eqbm_allaware]
    yvals_2 = [res[ks_2[1]]; [(res[ks_2[k]]*sum(shares[ks_2[k]][rows,:])-res[ks_2[k-1]]*sum(shares[ks_2[k-1]][rows,:]))/sum(shares[ks_2[k]][rows,:]-shares[ks_2[k-1]][rows,:]) for k=2:length(ks_2)]; ]
    #
    ks = [:baseline_by_zg_80,:direct_by_zg_80,:info_by_zg_80,:by_zg_80]
    yvals_3 = [res[ks[1]]; [(res[ks[k]]*sum(shares[ks[k]][rows,:])-res[ks[k-1]]*sum(shares[ks[k-1]][rows,:]))/sum(shares[ks[k]][rows,:]-shares[ks[k-1]][rows,:]) for k=2:length(ks)];]
    #
    yvals = [yvals_1; yvals_2; yvals_3]
    ll = min(0,min( minimum( yvals_1), minimum( yvals_2)) * 1.5)
    ul = min(5, max( maximum( yvals_1), maximum(yvals_2 )) * 1.75)
    pos = fig_2[nn,1]
    group = [fill(1,length(yvals_1)); 2; 3; fill(2,length(yvals_2)-2); fill(3,length(yvals_3))]
    ixx = [1; 3:5]
    ixx_aid = [1;2;2:5;]
    x = [ixx; ixx_aid; ixx;]
    xlab = [xlab3[ixx]; xlab3[ixx_aid]; xlab3[ixx];]
    mycolors = colors[group]; mycolors[6] = colors[2]
    barplot(pos, x, yvals ,axis = (xticks = (x, xlab), title="$name", 
    limits = (nothing,nothing,ll,ul),
    ),
    color_over_bar=:black,
    color_over_background=:black,
    bar_labels=[(y>0 ? string( round(y,digits=3)) : 
                 y<0 ? string( round(y,digits=3)) :
                 "") for y in yvals],
    #flip_labels_at = (minimum(yvals)+1, maximum(yvals)-1),
    dodge=group, color=mycolors, stack=x )
end
Makie.save(tablespath*"/fig_cf_1b_bar.png",fig_2)
#THIS SHOULD BE A TABLE? Or separate barcharts? 

# #complementary effects of information about aid and admissions
# xlab2 = ["baseline","full awareness","no TTT: eqbm.","TTT: direct","TTT: info","TTT: eqbm."]
# plots2 = []
# for (res,name) in zip([enroll_top10,enroll_rest,enroll_urm,enroll_theil,gpa_flagship,persist_flagship,stem_flagship,E_pi_flagship],
#     ["Top Decile (pct.)", "Not Top Decile (pct.)", "URM (pct.)", "Theil Index", "GPA", "1(Persist)", "1(STEM)","E(pi)"])
#     ys = [res[:baseline], res[:nottt_allaware], res[:base_allaware], res[:direct_allaware], res[:info_allaware], res[:eqbm_allaware]]
#     plt = waterfall(1:5,[ ys[2:end].-ys[1:end-1]; ],axis = (xticks = (1:5, xlab2[2:end]), title="Flagship Enrollment: $name"),)
#     push!(plots2,plt)
# end
keys_baseline = collect( keys(sumstats[:enroll_top10]))[[occursin("baseline",k) for k in String.(keys(sumstats[:enroll_top10]))]]
sumstats[:shareinfo] = OrderedDict( k=> 
(sumstats[:enroll_top10][:info] - sumstats[:enroll_top10][:direct]) /
(sumstats[:enroll_top10][:info] - sumstats[:enroll_top10][k]) for k in keys_baseline)

#limits of info (1): vary share autoadmitted
xlab3 = ["Baseline(0)","(5)", "(10)","(15)","(20)"]
fig_cf_3 = Figure(size=(1000,800))
for (nn,_res,name) in zip(1:8,[:enroll_top10,:enroll_rest,:enroll_urm,:enroll_theil,:gpa_flagship,:persist_flagship,:stem_flagship,:E_pi_flagship],
    ["Flagship Enrollment: Top Decile (% Change)", "Flagship Enrollment: Not Top Decile (% Change)", "Flagship Enrollment: URM (% Change)", "Flagship Enrollment: HS Concentration (Theil Index)", "GPA", "1(Persist)", "1(STEM)", "E(pi)"])
    res = sumstats[_res]
    pos = fig_cf_3[div(nn+1,2),mod(nn+1,2)+1]
    ys = [res[:baseline], res[:by_zg_95], res[:by_zg_90], res[:by_zg_85], res[:by_zg_80],]
    if nn <= 4
        yvals = 100* (ys[2:end].-ys[1:end-1])/ys[1]
    elseif nn==6 || nn==7
        yvals = 100* (ys[2:end].-ys[1:end-1])
    else
        yvals = (ys[2:end].-ys[1:end-1])
    end
    ll = min( 0, minimum(cumsum(yvals))) .* 1.5
    ul = max(0, maximum(cumsum(yvals))) .* 1.5
    _ax,_wf = waterfall(pos,1:4,yvals,axis = (xticks = (1:4, xlab3[2:end]), title="$name", limits=(nothing,nothing,ll,ul)), 
    #  bar_labels=[(y>0 ? "+"*string( round(y,digits=2))*"%" : 
    #  y<0 ? string( round(y,digits=2))*"%" :
    #  "") for y in yvals],
     )
    #hack to get waterfall plot to display bar labels
    bar_labels=[(y>0 ? "+"*string( round(y,digits=2))*"%" : 
    y<0 ? string( round(y,digits=2))*"%" :
    "") for y in yvals]
    wf_args = _wf.plots[1].args
    wf_attr = _wf.plots[1].attributes
    delete!(_ax,_wf)
    barplot!(_ax,wf_args...; wf_attr..., bar_labels=bar_labels,
    color_over_background=:black, color_over_bar=:black,)
end
Makie.save(tablespath*"/fig_cf_3.png",fig_cf_3)

#########################################
#table: direct and informative effects
#########################################
supertitles = "& \\multicolumn{4}{r}{{Pr(Enroll in Flagship), Treated Group}} & \\multicolumn{4}{c}{{Average Effect on Outcomes}} \\\\"

function print_df(df_point,df_se,titles; supertitles=supertitles, f=stdout, nleft=1, classbalance=false)
    rvec = repeat("c",length(titles)-nleft)
    lvec = repeat("l",nleft)
    println(f,"\\begin{tabular}{$(lvec)$(rvec)}")
    supertitles !== nothing && println(f, supertitles)
    println(f, join(titles, " & ") * " \\\\")
    println(f,"\\hline")
    start_at = classbalance ? 7 : 1
    for r = start_at:nrow(df_point)
        r==1 && println(f,"\\multicolumn{$(length(titles))}{l}{\\emph{A. Texas Top Ten}} \\\\ \\noalign{\\vskip 1mm}")
        r==3 && println(f,"\\multicolumn{$(length(titles))}{l}{\\emph{B. Automatic Admission By \$ z^\\text{admit}\\gamma \$}} \\\\ \\noalign{\\vskip 1mm}")
        r< 7 || classbalance || break
        r==7 && println(f,"\\multicolumn{$(length(titles))}{l}{\\emph{Class Balance: Additional Weight on Acad. Index}} \\\\ \\noalign{\\vskip 1mm}")
        println(f, join(string.(myround.(Vector( df_point[r,:]), digits=2)), " & ") * " \\\\")
        println(f, " & ("*join(myround.(Vector( df_se[r,2:end]),digits=2), ") & (") * ") \\\\ \\noalign{\\vskip 2mm} ")
    end
    println(f, "\\hline \\hline")
    println(f,"\\end{tabular}")
end

dfr = make_robustness_table()
dfr_point = cleanup_dfr(dfr)
dfr_point.cf[2] = "(2) +Aware"

dfr_se = deepcopy(dfr_point)
for nn in names(dfr_point)[2:end]
    dfr_se[!,Symbol(nn)] = std([dfrb[!,Symbol(nn)] for dfrb in dfr_boot])
end
dfr_se = cleanup_dfr(dfr_se)

dfr_titles = ["", "Base", "+ Mech.", "+ Info", "\\% Info", "\$\\Delta\$GPA",
    "\$\\Delta\$Persist","\$\\Delta\$STEM","\$\\Delta\$Payoff"]
print_df(dfr_point,dfr_se,dfr_titles; f=stdout, nleft=1, classbalance=false)

print_df(dfr_point,dfr_se,dfr_titles; f=stdout, nleft=1, classbalance=true)


open(tablespath*"/tab_robustness.tex","w") do f
    print_df(dfr_point,dfr_se,dfr_titles; f=f, nleft=1, classbalance=false)
end

open(tablespath*"/tab_robustness_classbalance_extension.tex","w") do f
    supertitles = "& \\multicolumn{4}{r}{{Pr(Enroll in Flagship), Treated Group}} & \\multicolumn{3}{c}{{Average Effect on Outcomes}} \\\\"
    print_df(dfr_point[:,1:end-1],dfr_se[:,1:end-1],dfr_titles[1:end-1]; f=f, nleft=1, supertitles=supertitles, classbalance=true)
end


#sanity check
meas_moved =  (sumstats[:enroll_top10][:eqbm] - sumstats[:enroll_top10][:baseline])/100  #2% of total pop, 0.02 
100 * (sumstats[:stem_flagship][:eqbm] - sumstats[:stem_flagship][:baseline])/meas_moved



#############################
#additional robustness table
#############################


dfr2a = make_robustness_table(outcomevar=:enroll_urm, academic_outcomes=false, _sum_wt=sum(popwt[urm]))
dfr2b = make_robustness_table(outcomevar=:enroll_highpoverty, academic_outcomes=false, _sum_wt=sum(popwt[highpoverty]))
dfr2c = make_robustness_table(outcomevar=:enroll_affluent, academic_outcomes=false, _sum_wt=sum(popwt[affluent]))
dfr2d = make_robustness_table(outcomevar=:enroll_theil, academic_outcomes=false, _sum_wt=1)
dfr_extra = [dfr2a,dfr2b,dfr2c,dfr2d]

colors = Makie.wong_colors()
fig1_alt = Figure(size=(1000,800))
mytitles = ["TTP", "TTP + Aware", "Top 20% zg", "4x Acad. Index"]
names_fig1_alt = [
    "A. Flagship Enrollment %: URM (Baseline=$(basevals[2])%)",
    "B. Flagship Enrollment %: High-Poverty High Schools (Baseline=$(basevals[3])%)",
    "C. Flagship Enrollment %: Affluent High Schools (Baseline=$(basevals[4])%)",
    "D. Flagship Enrollment: High School Concentration (Theil Index, Baseline=$(basevals[5]))"]

for rr=1:4
    (ll,ul) = (rr==1 ? (-0.7,3.1) : rr==2 ? (-1.1,3.8) : rr==3 ? (-3.5, 2.5) : (-0.1,0.1) )
    Label(fig1_alt[2rr-1,1:4], text=names_fig1_alt[rr], fontsize=16)
    for (cc,ctab) in enumerate([1,2,6,7])
        if cc == 2
            myrow = dfr_extra[rr][ctab,:]
            xticks = (1:5, ["Aid1","Aid2","Mech.","Info.","Eq."])
            myvec = [myrow[:info_aid] - dfr_extra[rr][1,:base], myrow[:base]-myrow[:info_aid], myrow[:direct], myrow[:info], myrow[:eqbm]]
            pos = Axis(fig1_alt[2rr,cc], xticks = xticks, title=mytitles[cc], limits=(nothing,nothing,ll,ul), yticklabelsvisible = false)
            waterfall!(pos,1:5, myvec, show_direction=true, show_final=true,color=colors[[4,5,1,2,3]],)
        else
            xticks = (1:3, ["Mech.","Info.","Eq."])
            myvec = Array(dfr_extra[rr][ctab,[:direct,:info,:eqbm]])
            pos = Axis(fig1_alt[2rr,cc], xticks = xticks, title=mytitles[cc], limits=(nothing,nothing,ll,ul), yticklabelsvisible = (cc==1))
            waterfall!(pos,1:3, myvec, show_direction=true, show_final=true,color=colors[1:3],)
        end 
    end
end
Legend(fig1_alt[9,:], [PolyElement(polycolor=colors[i]) for i=1:5], ["Mechanical Effect", "Informative Effect", "Equilibrium Effect", "Aid Info: Direct", "Aid Info: Eqbm."], orientation=:horizontal)
Makie.save(tablespath*"/fig_cf_1_alt.png",fig1_alt)









#############################
#numbers for paper 
#############################

# pct impact, "direct"
enroll_top10 = sumstats[:enroll_top10]
enroll_top10[:baseline] #baseline
(enroll_top10[:eqbm]-enroll_top10[:baseline]) #tot impact 
(enroll_top10[:eqbm]-enroll_top10[:baseline])/enroll_top10[:baseline] #pct
(enroll_top10[:info]-enroll_top10[:direct])/(enroll_top10[:eqbm] - enroll_top10[:baseline]) #what pct due to info
(enroll_top10[:direct]-enroll_top10[:baseline])/(enroll_top10[:eqbm] - enroll_top10[:baseline])

enroll_theil = sumstats[:enroll_theil]
enroll_theil[:baseline] #baseline
(enroll_theil[:eqbm]-enroll_theil[:baseline]) #tot impact 
(enroll_theil[:eqbm]-enroll_theil[:baseline])/enroll_theil[:baseline] #pct
(enroll_theil[:info]-enroll_theil[:direct])/(enroll_theil[:eqbm] - enroll_theil[:baseline])
(enroll_theil[:direct]-enroll_theil[:baseline])/(enroll_theil[:eqbm] - enroll_theil[:baseline])

