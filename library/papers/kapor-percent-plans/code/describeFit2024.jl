
################################################################
# 6. fit by class rank
################################################################
my_stats = getPrABC(zeros(J),1:J)
classrank = ziRaw[:,:classrank]
my_stats_killTTT = getPrABC(zeros(J),1:J, CC, classrank.==0.2, classrank.==0.2, θstar, popwt, 1.0, 1:II)

if use_old_simdata_approach
    #in this case, these files have not been loaded because admindata_new.jl was not run; load them now
    tamu = CSV.read(datapath*"/tamu.csv",DataFrame)
    uta = CSV.read(datapath*"/uta.csv",DataFrame)
    uta_rdsample = uta[ (uta.Ci .== 1), :] #drop missing classrank, classrank in bottom 10% of class
    tamu_rdsample = tamu[ (tamu.Ci .== 1), :] #drop missing classrank, classrank in bottom 10% of class
end

tamu_app = tamu[(Missings.replace(tamu[!,:classrank],1.0) .< .9), :]
uta_app = uta[(Missings.replace(uta[!,:classrank],1.0) .< .9), :]

gpa_hat_UTA = zq_UTA*dfoutcomes.gpa_uta
persist_hat_UTA = zq_UTA*dfoutcomes.persist_uta
stem_hat_UTA = zq_UTA*dfoutcomes.stem_uta
gpa_hat_TAMU = zq_TAMU*dfoutcomes.gpa_tamu
persist_hat_TAMU = zq_TAMU*dfoutcomes.persist_tamu
stem_hat_TAMU = zq_TAMU*dfoutcomes.stem_tamu


function summarize_mystats(my_stats,popwt=popwt,popwt_wv2=popwt_wv2,newinds=1:II,classrank=classrank[newinds])
    #_cr = sort( unique(classrank), rev=true)
    _cr = 0.9 : -0.1 : 0
    classr = round.(classrank,digits=1)
    myinds = [classr .== _c for _c in _cr]
    # use_old_simdata_approach || @assert sum(sum(myinds)) .== II   #in real data we have to get every obs; in simulations we have some GPA==1.0 which we can ignore here.
    model = OrderedDict(
        :prA => [[sum( my_stats[:A][j,inds].*popwt[inds]) ./ sum(popwt[inds]) for inds in myinds] for j=1:J],
        :prB => [[sum(my_stats[:B][j,inds].*popwt[inds]) ./ sum(popwt[inds]) for inds in myinds] for j=1:J],
        :prB!A => [[sum(my_stats[:B][j,inds].*popwt[inds]) ./ sum(my_stats[:A][j,inds].*popwt[inds]) for inds in myinds] for j=1:J],
        :prC!B => [[sum(my_stats[:C][j,inds].*popwt[inds]) ./ sum(my_stats[:B][j,inds].*popwt[inds]) for inds in myinds] for j=1:J],
        :prC => [[sum(my_stats[:C][j,inds].*popwt[inds])./sum(popwt[inds]) for inds in myinds] for j=1:J],
        :prPersist!C => let 
            outc = Dict(6=>persist_hat_TAMU, 7=>persist_hat_UTA)
            [ [sum(my_stats[:C][j,inds] .* outc[j][inds].*popwt[inds]) ./ sum(my_stats[:C][j,inds].*popwt[inds]) for inds in myinds] for j=6:J] 
        end,
        :prSTEM!C => let 
            outc = Dict(6=>stem_hat_TAMU, 7=>stem_hat_UTA)
            [ [sum(my_stats[:C][j,inds] .* outc[j][inds].*popwt[inds]) ./ sum(my_stats[:C][j,inds].*popwt[inds]) for inds in myinds] for j=6:J] 
        end,   
        :GPA!C => let 
            outc = Dict(6=>gpa_hat_TAMU, 7=>gpa_hat_UTA)
            [ [sum(my_stats[:C][j,inds] .* outc[j][inds].*popwt[inds]) ./ sum(my_stats[:C][j,inds].*popwt[inds]) for inds in myinds] for j=6:J] 
        end,
    )
    for k in [:prPersist!C,:prSTEM!C, :GPA!C]
        mk = model[k]
        push!(mk, (mk[1].*model[:prC][6] .+ mk[2].*model[:prC][7])./(model[:prC][6] + model[:prC][7]))
    end
    model, myinds
end
    
model, myinds = summarize_mystats(my_stats)
model_killTTT, _ = summarize_mystats(my_stats_killTTT)

prA_surv = [let
    tmp = [CC.portf[j,CC.Aij[ii]] .> 0 for ii=1:II] 
    [mean(tmp[inds], Weights(popwt[inds])) for inds in myinds]
end for j=1:J]
prB_surv = [let
    prb = [CC.portf[j,CC.Bij[ii]] .> 0 for ii=1:II]
    [mean(prb[inds],Weights(popwt[inds])) for inds in myinds]
end for j=1:J]
prB!A_surv = [let
    pra = [CC.portf[j,CC.Aij[ii]] .> 0 for ii=1:II] .* popwt
    prb = [CC.portf[j,CC.Bij[ii]] .> 0 for ii=1:II] .* popwt
    [sum(prb[inds])/sum(pra[inds]) for inds in myinds]
end for j=1:J]
prC!B_surv = [let
    prc = (CC.Ci.==j) .* popwt_wv2
    prb = [CC.portf[j,CC.Bij[ii]] .> 0 for ii=1:II] .* popwt_wv2
    [sum(prc[inds])/sum(prb[inds]) for inds in myinds]
end for j=1:J]

fullappsample = [uta_app; tamu_app]
fullrdsample = [tamu_rdsample; uta_rdsample]
prB!A_admin = [[mean(df.Bij[lb .< df.classrank .<= lb+.1]) for lb = 0.9: -0.1 :0] for df in [tamu_app, uta_app, fullappsample]]
prC!B_admin = [[mean(df.Ci[(df.Bij.==1) .& (lb .< df.classrank .<= lb+.1)]) for lb = 0.9: -0.1 :0] for df in [tamu_app, uta_app, fullappsample]]
prB_admin = [prB!A_admin[1] .* prA_surv[6], prB!A_admin[2] .* prA_surv[7], prB!A_admin[1] .* prA_surv[6] .+ prB!A_admin[2] .* prA_surv[7]]
prC_admin = [prB_admin[j] .* prC!B_admin[j] for j=1:3]
prPersist!C_admin = [[mean(df.persist[lb .< df.classrank .<= lb+.1]) for lb = 0.9: -0.1 :0] for df in [tamu_rdsample, uta_rdsample, fullrdsample]]
prSTEM!C_admin = [[mean(df.stem[lb .< df.classrank .<= lb+.1]) for lb = 0.9: -0.1 :0] for df in [tamu_rdsample, uta_rdsample, fullrdsample]]
GPA!C_admin = [[mean(df.gpa[lb .< df.classrank .<= lb+.1]) for lb = 0.9: -0.1 :0] for df in [tamu_rdsample, uta_rdsample, fullrdsample]]

cellcount_apply_admin = [[sum(lb .< df.classrank .<= lb+.1) for lb = 0.9: -0.1 :0] for df in [tamu_app, uta_app, fullappsample]]
cellcount_admit_admin = [[sum( df.Bij .* (lb .< df.classrank .<= lb+.1)) for lb = 0.9: -0.1 :0] for df in [tamu_app, uta_app, fullappsample]]
cellcount_enroll_admin = [[sum(lb .< df.classrank .<= lb+.1) for lb = 0.9: -0.1 :0] for df in [tamu_rdsample, uta_rdsample, fullrdsample]]



#panel 1: pr(apply UTA), model vs. survey. also do apply TAMU; apply any public; apply any 4-year
#panel 2: pr(admit | apply UTA) by class rank, model vs. survey vs. admin data
#panel 3: pr(enroll UTA | admit UTA) by class rank, model v. survey v. admin. also do pr(enroll TAMU;); pr(enroll any public); pr(enroll anywhere)
#panel 4: pr(persist UTA | enroll UTA) by class rank, model vs. admin
#panel 5: pr(stem UTA | enroll UTA) by class rank, model vs. admin
#panel 6: pr(gpa UTA | enroll UTA) by class rank, model vs. admin
ix_surv = 6:7 

function make_fitfig(ix_surv=6:7)
    min_n = 20
    ix_admin =  ix_surv==(6:7) ? 3 : (ix_surv==(6:6)) ? 1 : (ix_surv==(7:7)) ? 2 : 0
    _ncs = [c in ix_surv for c in CC.Ci] #(CC.Ci.==7)
    _nbs = [any( CC.portf[ix_surv,CC.Bij[ii]] .> 0) for ii=1:II]
    _nas = [any( CC.portf[ix_surv,CC.Aij[ii]] .> 0) for ii=1:II]
    vcs = [sum(_ncs[inds]) for inds in myinds] .> min_n #survey cell counts
    vbs = [sum(_nbs[inds]) for inds in myinds] .> min_n
    vas = [sum(_nas[inds]) for inds in myinds] .> min_n
    vaa = cellcount_apply_admin[ix_admin] .> min_n
    vba = cellcount_admit_admin[ix_admin] .> min_n
    vca = cellcount_enroll_admin[ix_admin] .> min_n

    fig_insample = Figure(size=(1000,700))
    rgrid = 0.05:0.1:.95
    mycolors = [:blue,:orange,:green,:purple]
    mytag = ix_surv==(6:7) ? "Flagship Universities" : ix_surv==(6:6) ? "Texas A&M" : "UT Austin"
    pos = Axis(fig_insample[1,1], xlabel = "Class Rank", title="A. Num. Apps, $mytag")
    scatter!(pos,rgrid,sum(model[:prA][ix_surv]), color = mycolors[1], marker=:circle)
    scatter!(pos,rgrid[vas],sum(prA_surv[ix_surv])[vas], color = mycolors[2], marker=:diamond)
    scatter!(pos,rgrid[end:end],sum(model_killTTT[:prA][ix_surv])[end:end], color = mycolors[4], marker=:xcross)

    pos = Axis(fig_insample[1,2], xlabel = "Class Rank", title="B. Num. Offers, $mytag")
    bmodel = sum([model[:prB!A][j] .* model[:prA][j] for j in ix_surv]) #./ sum(model[:prA][ix_surv])
    bsurv = sum([prB!A_surv[j] .* prA_surv[j] for j in ix_surv]) #./ sum(prA_surv[ix_surv])
    bmodel_killTTT = sum([model_killTTT[:prB!A][j] .* model_killTTT[:prA][j] for j in ix_surv]) #./ sum(model_killTTT[:prA][ix_surv])
    badmin = prB_admin[ix_admin] #prB!A_admin[ix_admin] 
    scatter!(pos,rgrid,bmodel, color = mycolors[1], marker=:circle)
    scatter!(pos,rgrid[vbs],bsurv[vbs], color = mycolors[2], marker=:diamond)
    scatter!(pos,rgrid[vba],badmin[vba], color = mycolors[3], marker=:rect)
    scatter!(pos,rgrid[end:end],bmodel_killTTT[end:end], color = mycolors[4], marker=:xcross)

    pos = Axis(fig_insample[2,1], xlabel = "Class Rank", title="C. Pr(Enroll), $mytag")
    cmodel = sum([model[:prC!B][j] .* model[:prB][j] for j in ix_surv]) #./ sum(model[:prB][ix_surv])
    csurv = sum([prC!B_surv[j] .* prB_surv[j] for j in ix_surv]) #./ sum(prB_surv[ix_surv])
    cmodel_killTTT = sum([model_killTTT[:prC!B][j] .* model_killTTT[:prB][j] for j in ix_surv]) #./ sum(model_killTTT[:prB][ix_surv])
    cadmin = prC_admin[ix_admin] #prC!B_admin[ix_admin]
    scatter!(pos,rgrid,cmodel, color = mycolors[1], marker=:circle)
    scatter!(pos,rgrid[vcs],csurv[vcs], color = mycolors[2], marker=:diamond)
    scatter!(pos,rgrid[vca],cadmin[vca], color = mycolors[3], marker=:rect)
    scatter!(pos,rgrid[end:end],cmodel_killTTT[end:end], color = mycolors[4], marker=:xcross)

    pos = Axis(fig_insample[2,2], xlabel = "Class Rank", title="D. Pr(Persist | Enroll), $mytag")
    scatter!(pos,rgrid,model[:prPersist!C][ix_admin], color = mycolors[1], marker=:circle)
    scatter!(pos,rgrid[vca],prPersist!C_admin[ix_admin][vca], color = mycolors[3], marker=:rect)
    scatter!(pos,rgrid[end:end],model_killTTT[:prPersist!C][ix_admin][end:end], color = mycolors[4], marker=:xcross)

    pos = Axis(fig_insample[3,1], xlabel = "Class Rank", title="E. Pr(STEM | Enroll), $mytag")
    scatter!(pos,rgrid,model[:prSTEM!C][ix_admin], color = mycolors[1], marker=:circle)
    scatter!(pos,rgrid[vca],prSTEM!C_admin[ix_admin][vca], color = mycolors[3], marker=:rect)
    scatter!(pos,rgrid[end:end],model_killTTT[:prSTEM!C][ix_admin][end:end], color = mycolors[4], marker=:xcross)

    pos = Axis(fig_insample[3,2], xlabel = "Class Rank", title="F. GPA | Enroll, $mytag")
    scatter!(pos,rgrid,model[:GPA!C][ix_admin], color = mycolors[1], marker=:circle)
    scatter!(pos,rgrid[vca],GPA!C_admin[ix_admin][vca], color = mycolors[3], marker=:rect)
    scatter!(pos,rgrid[end:end],model_killTTT[:GPA!C][ix_admin][end:end], color = mycolors[4], marker=:xcross)

    mymarkers = [:circle,:diamond,:rect,:xcross]
    #el_fit = [PolyElement(polycolor=colors[i]) for i=1:length(mycolors)]
    el_fit = [MarkerElement(color=colors[i], marker=mymarkers[i]) for i=1:length(mycolors)]
    Legend(fig_insample[4,1:2], el_fit, ["Model","Survey","Admin. Data","Model: Drop TTT"], orientation=:horizontal)
    fig_insample
end

function make_fitfig_survOnly()
    min_n = 15
    fig_insample = Figure(size=(1000,700))
    rgrid = 0.05:0.1:.95
    mycolors = [:blue,:orange,:green,:purple]
    for (m,ix_surv,mytag) in zip(1:2,[1:7,findall(JTXpublic)],["All 4-year", "All Public"])
        _ncs = [c in ix_surv for c in CC.Ci] #(CC.Ci.==7)
        _nbs = [any( CC.portf[ix_surv,CC.Bij[ii]] .> 0) for ii=1:II]
        _nas = [any( CC.portf[ix_surv,CC.Aij[ii]] .> 0) for ii=1:II]
        vcs = [sum(_ncs[inds]) for inds in myinds] .> min_n #survey cell counts
        vbs = [sum(_nbs[inds]) for inds in myinds] .> min_n
        vas = [sum(_nas[inds]) for inds in myinds] .> min_n
        pos = Axis(fig_insample[1,m], xlabel = "Class Rank", title="A. Num. Apps, $mytag")
        scatter!(pos,rgrid,sum(model[:prA][ix_surv]), color = mycolors[1], marker=:circle)
        scatter!(pos,rgrid[vas],sum(prA_surv[ix_surv])[vas], color = mycolors[2], marker=:diamond)
        #scatter!(pos,rgrid[end:end],sum(model_killTTT[:prA][ix_surv])[end:end], color = mycolors[4], marker=:xcross)
        pos = Axis(fig_insample[2,m], xlabel = "Class Rank", title="B. Num. Offers, $mytag")
        bmodel = sum([model[:prB!A][j] .* model[:prA][j] for j in ix_surv]) #./ sum(model[:prA][ix_surv])
        bsurv = sum([prB!A_surv[j] .* prA_surv[j] for j in ix_surv]) #./ sum(prA_surv[ix_surv])
        bmodel_killTTT = sum([model_killTTT[:prB!A][j] .* model_killTTT[:prA][j] for j in ix_surv]) ./ sum(model_killTTT[:prA][ix_surv])
        scatter!(pos,rgrid,bmodel, color = mycolors[1], marker=:circle)
        scatter!(pos,rgrid[vbs],bsurv[vbs], color = mycolors[2], marker=:diamond)
        #scatter!(pos,rgrid[end:end],bmodel_killTTT[end:end], color = mycolors[4], marker=:xcross)
        pos = Axis(fig_insample[3,m], xlabel = "Class Rank", title="C. Pr(Enroll), $mytag")
        cmodel = sum([model[:prC!B][j] .* model[:prB][j] for j in ix_surv]) #./ sum(model[:prB][ix_surv])
        csurv = sum([prC!B_surv[j] .* prB_surv[j] for j in ix_surv]) #./ sum(prB_surv[ix_surv])
        cmodel_killTTT = sum([model_killTTT[:prC!B][j] .* model_killTTT[:prB][j] for j in ix_surv]) ./ sum(model_killTTT[:prB][ix_surv])
        scatter!(pos,rgrid,cmodel, color = mycolors[1], marker=:circle)
        scatter!(pos,rgrid[vcs],csurv[vcs], color = mycolors[2], marker=:diamond)
        #scatter!(pos,rgrid[end:end],cmodel_killTTT[end:end], color = mycolors[4], marker=:xcross)
    end
    mymarkers = [:circle,:diamond]
    #el_fit = [PolyElement(polycolor=colors[i]) for i=1:length(mycolors)]
    el_fit = [MarkerElement(color=colors[i], marker=mymarkers[i]) for i=1:length(mymarkers)]
    Legend(fig_insample[4,1:2], el_fit, ["Model","Survey"], orientation=:horizontal)
    fig_insample
end

fig_insample = make_fitfig(6:7)
Makie.save(tablespath*"/fig_fit_1.png",fig_insample)

fig_insample_UTA = make_fitfig(7:7)
fig_insample_TAMU = make_fitfig(6:6)
Makie.save(tablespath*"/fig_fit_1_UTA.png",fig_insample_UTA)
Makie.save(tablespath*"/fig_fit_1_TAMU.png",fig_insample_TAMU)

fig_insample_other = make_fitfig_survOnly()
Makie.save(tablespath*"/fig_fit_1_Other.png",fig_insample_other)















################################################################
# 7. fit - table
################################################################
#old table w/ total effects; show top-decile and urm flagship enrollment shares in 1997, 2001, 2002.

################################################################
# 8. numbers for paper
################################################################
pr_flagship_topdec = model[:prC][7][end] + model[:prC][6][end]
pr_flagship_topdec_killTTT = model_killTTT[:prC][7][end] + model_killTTT[:prC][6][end]
pr_flagship_dif = pr_flagship_topdec - pr_flagship_topdec_killTTT

open(numberspath*"/fit_topdec_share_model.tex","w") do f
    print(f, round(100*pr_flagship_topdec, digits=1))
end
open(numberspath*"/fit_rdd.tex","w") do f
    print(f, round(100*pr_flagship_dif, digits=1))
end

jump_prA_model_UTA = model[:prA][7][end] - model[:prA][7][end-1]
jump_prA_surv_UTA = prA_surv[7][end] - prA_surv[7][end-1]

jump_prA_model_TAMU = model[:prA][6][end] - model[:prA][6][end-1]
jump_prA_surv_TAMU = prA_surv[6][end] - prA_surv[6][end-1]

jump_prA_model = jump_prA_model_UTA + jump_prA_model_TAMU
jump_prA_surv = jump_prA_surv_UTA + jump_prA_surv_TAMU

open(numberspath*"/jump_prA_surv.tex","w") do f
    print(f, round(100*jump_prA_surv, digits=1))
end
open(numberspath*"/jump_prA_model.tex","w") do f
    print(f, round(100*jump_prA_model, digits=1))
end
open(numberspath*"/prA_model.tex","w") do f
    print(f, round(100*(model[:prA][7][end] + model[:prA][6][end]), digits=1))
end

open(numberspath*"/prA_model_dif_nottt.tex","w") do f
    print(f, round(100*(model[:prA][7][end] - model_killTTT[:prA][7][end] + model[:prA][6][end] - model_killTTT[:prA][6][end]), digits=1))
end

open(numberspath*"/jump_prA_surv_UTA.tex","w") do f
    print(f, round(100*jump_prA_surv_UTA, digits=1))
end
open(numberspath*"/jump_prA_model_UTA.tex","w") do f
    print(f, round(100*jump_prA_model_UTA, digits=1))
end
open(numberspath*"/prA_model_UTA.tex","w") do f
    print(f, round(100*model[:prA][7][end], digits=1))
end

open(numberspath*"/prA_model_dif_nottt_UTA.tex","w") do f
    print(f, round(100*(model[:prA][7][end] - model_killTTT[:prA][7][end]), digits=1))
end

