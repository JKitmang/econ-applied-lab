##################################
# code to reuse gradient temp storage
##################################

function myobj_nottt(Δπ,tot_enroll_ttt=tot_enroll_ttt,CC=CC,θ=θstar,ttta=falses(II),tttb=falses(II),γq=1.0,IX=1:II,popwt=popwt)
    tot_enroll_nottt = sum( getPrEnroll(Δπ, findall(JTXpublic), CC, ttta, tttb, θ, popwt, γq, IX))
    dif = (tot_enroll_ttt - tot_enroll_nottt).*10000
    println("Δπ = $(ForwardDiff.value.(Δπ)), dif=$(ForwardDiff.value.(dif))")
    val = sum( dif.^2)
end

#callable struct in order to reuse type for gradient
struct CF{T}
    cap::Array{Float64,1}
    CC::T
    ttta::BitVector
    tttb::BitVector
    θ::Array{Float64,1}
    γq::Float64
    inds::Array{Int,1}
    popwt::Array{Float64,1}
end
CF(cap,CC::T,ttta,tttb,θ,γq) where {T} = CF{T}(cap,CC,ttta,tttb,θ,γq,1:II,popwt)

struct DCF{T}
    cf::CF{T}
end
(cf::CF{T})(Δπ) where {T} = myobj_nottt(Δπ,cf.cap,cf.CC,cf.θ,cf.ttta,cf.tttb,cf.γq,cf.inds,cf.popwt)
(dcf::DCF{T})(dst,Δπ) where {T} = (dst .= ForwardDiff.gradient(dcf.cf,Δπ,ForwardDiff.GradientConfig(nothing,zeros(3),ForwardDiff.Chunk{3}())))

##################################
##################################

function getPrABC(Δπ_inds, inds_C=1:J, CC=CC, ttti_app=CC.ttt, ttti_admit=CC.ttt, θ=θstar, popwt=popwt, γq=1.0, IX=1:II)
    T = eltype(Δπ_inds)
    Δπ = zeros(T,J)
    Δπ[inds_C] .= Δπ_inds
    II = length(IX)
    dsts = OrderedDict(:A=>zeros(7,II), :B=>zeros(7,II), :C=>zeros(7,II), :A_public => zeros(II), :B_public=>zeros(II), 
        :A_all=>zeros(II), :B_all=>zeros(II), :A_flagship=>zeros(II), :B_flagship=>zeros(II))
    Threads.@threads for ii=1:II
        other_stats = Dict(:apply=>zeros(T,7), :admit=>zeros(T,7), :apply_admit_summary => zeros(T,3,2))
        dsts[:C][:,ii] .= getPrEnroll_i(Δπ, ii, inds_C, CC, ttti_app[ii], ttti_admit[ii], θ, γq, IX[ii], other_stats=other_stats) 
        dsts[:B][:,ii] .= other_stats[:admit] 
        dsts[:A][:,ii] .= other_stats[:apply] 
        dsts[:A_public][ii] = other_stats[:apply_admit_summary][1,1]
        dsts[:A_all][ii] = other_stats[:apply_admit_summary][2,1]
        dsts[:A_flagship][ii] = other_stats[:apply_admit_summary][3,1]
        dsts[:B_public][ii] = other_stats[:apply_admit_summary][1,2]
        dsts[:B_all][ii] = other_stats[:apply_admit_summary][2,2]
        dsts[:B_flagship][ii] = other_stats[:apply_admit_summary][3,2]
    end
    dsts
end

function getPrEnroll(Δπ_inds, inds_C=1:J, CC=CC, ttti_app=CC.ttt, ttti_admit=CC.ttt, θ=θstar, popwt=popwt, γq=1.0, IX=1:II)
    T = eltype(Δπ_inds)
    Δπ = zeros(T,J)
    Δπ[inds_C] .= Δπ_inds
    II = length(IX)
    dsts = [zeros(T,7) for ii=1:II]
    Threads.@threads for ii=1:II
        dsts[ii] .= getPrEnroll_i(Δπ, ii, inds_C, CC, ttti_app[ii], ttti_admit[ii], θ, γq, IX[ii]) .* popwt[ii]
    end
    dsts
end

function getPrEnroll_i(Δπ, ii, inds_C=1:J, CC=CC, ttti_app=CC.ttt[ii], ttti_admit=CC.ttt[ii], θ=θstar, γq=1.0, ix=ii; other_stats=Dict())
    T = eltype(Δπ)
    ws::WS{T} = get!(WSS[ii],T) do 
        WS(CC,ii,T)
    end
    #unpack parameters
    γadmit = view(θ,1:J+nz)  #first j elements are cutoffs; next nz are zi 
    γs = SVector{3}(θ[J+nz+1],θ[J+nz+2],θ[J+nz+3])
    γq!s = SVector{3}(θ[J+nz+4],θ[J+nz+5],θ[J+nz+6])
    λ = 1/(1+exp(-θ[J+nz+7])) #correlation parameter for final-period shocks
    γcfixed = SVector{3}(θ[J+nz+8],θ[J+nz+9],θ[J+nz+10]) #app costs
    γcvariable = SVector{3}(θ[J+nz+11],θ[J+nz+12],θ[J+nz+13])
    γAppshocks_0 = θ[J+nz+14]
    βp = SVector{2}(θ[J+nz+15],θ[J+nz+16])
    βx = view(θ,J+nz+17 : J+nz+16+nx) #deterministic coefficients
    log_σ_rc = view(θ,J+nz+17+nx: J+nz+16+nx+nrc)
    ws.σ_rc .= exp.(log_σ_rc)
    βx_aware = view(θ,J+nz+17+nx+nrc:J+nz+16+2nx+nrc)
    βy_aware = θ[J+nz+17+2nx+nrc]
    σ_βi0_aware = exp(θ[J+nz+18+2nx+nrc])
    coef_aidamount = view(θ,J+nz+18+2nx+nrc+1:2J+nz+18+2nx+nrc) #multiplies (Pj-efc_i); cost at j if aware.
    log_σ_e = view(θ,2J+nz+19+2nx+nrc : 3J+nz+18+2nx+nrc)
    γAppshocks_pov = θ[3J+nz+19+2nx+nrc]
    γAppshocks_urm = θ[3J+nz+20+2nx+nrc]
    σAppShocks = log(1 + exp(γAppshocks_0*CC.zi_info[ii,1] + γAppshocks_pov*CC.zi_info[ii,2] + γAppshocks_urm*CC.zi_info[ii,3]))
    #...
    #construct some useful terms
    var_s = log(1 + exp(dot(γs,view(CC.zi_info,ii,:))))
    var_q!s = log(1 + exp(dot(γq!s,view(CC.zi_info,ii,:))))
    # var_q = var_s + var_q!s
    # sig_q = sqrt(var_q)
    sig_q!s = sqrt(var_q!s)
    c0 = dot(γcfixed, view(CC.zi_info,ii,:))
    c1 = dot(γcvariable, view(CC.zi_info,ii,:))
    zγ = ws.zγ; mul!(zγ,CC.zi[ii],γadmit)  #all schools
    grid = CC.grid[ii]
    ###################
    # from here on is new
    ###################
    zγ .-= Δπ
    dst = zeros(T,J)
    for mm=1:ndraws
        # deal w/ likelihood of awareness draws
        lik_aware = getFinaidAwarenessLikelihood(view(grid.aware,:,mm),βx_aware,βy_aware,σ_βi0_aware,CC.xij[ii],grid.y[mm],grid.efc[mm],ws.temps_j,nothing)
        ws.ll[mm] = log(lik_aware) - log(grid.pr0_aware[mm])
        # first, get ex-ante app probabilities prA
        s = grid.s[mm]*sqrt(var_s)
        ui, aid_m = getu(βx,ws.σ_rc,log_σ_e,βp,coef_aidamount,CC,ii,mm,ws)
        mul!(ws.invT_UAd[mm],CC.invTT,getvaladmit(ui,λ,ws,CC))
        Pq = get_Pq(s,sig_q!s,zγ,ws,ii,CC,ttti_app,γq) #ws.Pq
        V = getvalapp(ws.invT_UAd[mm],Pq,ws,CC,mm) #also fills Ps
        V .-= getappcosts(c0,c1,CC,ii,ws) #deterministic part of costs
        prA = ws.ValApp_allM[mm]
        lik_logit_all!(prA,V,σAppShocks)
        # now integrate over offers and enrollment decisions
        if ttti_app == ttti_admit
            Ps = ws.Ps[mm] #filled by getvalapp; have already computed this object
        else #need to recompute Pq and Ps, probability of admission to all programs in each B
            Pq = get_Pq(s,sig_q!s,zγ,ws,ii,CC,ttti_admit,γq) #ws.Pq
            Ps = ws.Ps[mm] #probability of admission to sets of programs | s
            Ps .= 0
            for ind_q = 1:length(Qnodes)
                Ps .+= Pq[ind_q].*Qweights[ind_q]
            end
        end
        prA_TT_Ps_transpose = ws.temps_a[1]
        mul!(prA_TT_Ps_transpose,CC.TTp,prA) #mul!(prA_TT,prA',CC.TT)
        #
        prA_TT_Ps_transpose .*= Ps
        prA_TT_Ps_invTT_transpose = ws.temps_a[2]
        mul!(prA_TT_Ps_invTT_transpose,CC.invTT',prA_TT_Ps_transpose)
        # offers 
        for Ci in inds_C
            prj!B = ws.temps_a[1]
            prj!B .= 0
            for bb=1:nportf
                if CC.portf[Ci,bb] > 0
                    prj!B[bb] = getenrollmentprob(ui,λ,bb,Ci,CC,ws)
                end
            end
            dst[Ci] += dot(prA_TT_Ps_invTT_transpose,prj!B)*exp(ws.ll[mm])
        end
        if :apply in keys(other_stats)
            other_stats[:apply] .+= [dot(prA, (CC.portf[j,:].>0)) for j=1:J] .* exp(ws.ll[mm])
        end
        if :admit in keys(other_stats)
            prB = prA_TT_Ps_invTT_transpose
            other_stats[:admit] .+= [dot(prB, (CC.portf[j,:].>0)) for j=1:J] .* exp(ws.ll[mm])
        end
        if :apply_admit_summary in keys(other_stats)
            inds_public = vec(sum( CC.portf[JTXpublic,:], dims=1))
            inds_all = vec(sum(CC.portf,dims=1))
            inds_flagship = vec(maximum(CC.portf[6:7,:],dims=1))
            prB = prA_TT_Ps_invTT_transpose
            other_stats[:apply_admit_summary][1,1] += dot(prA,inds_public) * exp(ws.ll[mm])
            other_stats[:apply_admit_summary][2,1] += dot(prA,inds_all) * exp(ws.ll[mm])
            other_stats[:apply_admit_summary][3,1] += dot(prA,inds_flagship) * exp(ws.ll[mm])

            other_stats[:apply_admit_summary][1,2] += dot(prB,inds_public) * exp(ws.ll[mm])
            other_stats[:apply_admit_summary][2,2] += dot(prB,inds_all) * exp(ws.ll[mm])
            other_stats[:apply_admit_summary][3,2] += dot(prB,inds_flagship) * exp(ws.ll[mm])
        end
    end
    dst ./= sum(exp.(ws.ll))
    for k in keys(other_stats)
        other_stats[k] ./= sum(exp.(ws.ll))
    end
    dst
end


function lik_logit_all!(dst,v,σ)
    maxv = maximum(v)
    dst .= exp.( (v .- maxv)./σ)
    dst ./= sum(dst)
end




function getPosteriorQ(Δπ_inds, inds_C=1:J, CC=CC, ttti_app=CC.ttt, ttti_admit=CC.ttt, θ=θstar, popwt=popwt, γq=1.0, IX=1:II)
    T = eltype(Δπ_inds)
    Δπ = zeros(T,J)
    Δπ[inds_C] .= Δπ_inds
    II = length(IX)
    #dsts = [zeros(T,7) for ii=1:II]
    Eshock = zeros(T,J,II)
    Eq = zeros(T,J,II)
    Ev = zeros(T,J,II)
    Threads.@threads for ii=1:II
        #println("getPosteriorQ: $ii")
        Eqi,Eshocki,Evi = getPosteriorQ_i(Δπ, ii, inds_C, CC, ttti_app[ii], ttti_admit[ii], θ, γq ,IX[ii])
        Eq[:,ii] .= Eqi
        Eshock[:,ii] .= Eshocki
        Ev[:,ii] .= Evi
    end
    Eq,Eshock,Ev
end


function getPosteriorQ_i(Δπ, ii, inds_C=1:J, CC=CC, ttti_app=CC.ttt[ii], ttti_admit=CC.ttt[ii], θ=θstar, γq=1.0, ix=ii)
    T = eltype(Δπ)
    ws::WS{T} = get!(WSS[ii],T) do 
        WS(CC,ii,T)
    end
    #unpack parameters
    γadmit = view(θ,1:J+nz)  #first j elements are cutoffs; next nz are zi 
    γs = SVector{3}(θ[J+nz+1],θ[J+nz+2],θ[J+nz+3])
    γq!s = SVector{3}(θ[J+nz+4],θ[J+nz+5],θ[J+nz+6])
    λ = 1/(1+exp(-θ[J+nz+7])) #correlation parameter for final-period shocks
    γcfixed = SVector{3}(θ[J+nz+8],θ[J+nz+9],θ[J+nz+10]) #app costs
    γcvariable = SVector{3}(θ[J+nz+11],θ[J+nz+12],θ[J+nz+13])
    γAppshocks_0 = θ[J+nz+14]
    βp = SVector{2}(θ[J+nz+15],θ[J+nz+16])
    βx = view(θ,J+nz+17 : J+nz+16+nx) #deterministic coefficients
    log_σ_rc = view(θ,J+nz+17+nx: J+nz+16+nx+nrc)
    ws.σ_rc .= exp.(log_σ_rc)
    βx_aware = view(θ,J+nz+17+nx+nrc:J+nz+16+2nx+nrc)
    βy_aware = θ[J+nz+17+2nx+nrc]
    σ_βi0_aware = exp(θ[J+nz+18+2nx+nrc])
    coef_aidamount = view(θ,J+nz+18+2nx+nrc+1:2J+nz+18+2nx+nrc) #multiplies (Pj-efc_i); cost at j if aware.
    log_σ_e = view(θ,2J+nz+19+2nx+nrc : 3J+nz+18+2nx+nrc)
    γAppshocks_pov = θ[3J+nz+19+2nx+nrc]
    γAppshocks_urm = θ[3J+nz+20+2nx+nrc]
    σAppShocks = log(1 + exp(γAppshocks_0*CC.zi_info[ii,1] + γAppshocks_pov*CC.zi_info[ii,2] + γAppshocks_urm*CC.zi_info[ii,3]))
    #...
    #construct some useful terms
    var_s = log(1 + exp(dot(γs,view(CC.zi_info,ii,:))))
    var_q!s = log(1 + exp(dot(γq!s,view(CC.zi_info,ii,:))))
    # var_q = var_s + var_q!s
    # sig_q = sqrt(var_q)
    sig_q!s = sqrt(var_q!s)
    c0 = dot(γcfixed, view(CC.zi_info,ii,:))
    c1 = dot(γcvariable, view(CC.zi_info,ii,:))
    zγ = ws.zγ; mul!(zγ,CC.zi[ii],γadmit)  #all schools
    grid = CC.grid[ii]
    ###################
    # from here on is new
    ###################
    zγ .-= Δπ
    _denom_ = zeros(T,J)
    Eq = zeros(T,J)
    Eshock = zeros(T,J)
    Ev = zeros(T,J)
    _pr_ = zeros(T,J) #temp storage for enrollment probs
    for mm=1:ndraws
        # deal w/ likelihood of awareness draws
        lik_aware = getFinaidAwarenessLikelihood(view(grid.aware,:,mm),βx_aware,βy_aware,σ_βi0_aware,CC.xij[ii],grid.y[mm],grid.efc[mm],ws.temps_j,nothing)
        ws.ll[mm] = log(lik_aware) - log(grid.pr0_aware[mm])
        # first, get ex-ante app probabilities prA
        s = grid.s[mm]*sqrt(var_s)
        ui, aid_m = getu(βx,ws.σ_rc,log_σ_e,βp,coef_aidamount,CC,ii,mm,ws)
        mul!(ws.invT_UAd[mm],CC.invTT,getvaladmit(ui,λ,ws,CC))
        Pq = get_Pq(s,sig_q!s,zγ,ws,ii,CC,ttti_app,γq) #ws.Pq
        V = getvalapp(ws.invT_UAd[mm],Pq,ws,CC,mm) #also fills Ps
        V .-= getappcosts(c0,c1,CC,ii,ws) #deterministic part of costs
        prA = ws.ValApp_allM[mm]
        lik_logit_all!(prA,V,σAppShocks)
        # now integrate over offers and enrollment decisions
        Pq = get_Pq(s,sig_q!s,zγ,ws,ii,CC,ttti_admit,γq) #ws.Pq
        for (_qq,wq,Pqq) in zip(Qnodes,Qweights,Pq)
            qq = s + _qq*sqrt(2)*sig_q!s
            prA_TT_Ps_transpose = ws.temps_a[1]
            mul!(prA_TT_Ps_transpose,CC.TTp,prA) #mul!(prA_TT,prA',CC.TT)
            prA_TT_Ps_transpose .*= Pqq
            prA_TT_Ps_invTT_transpose = ws.temps_a[2]
            mul!(prA_TT_Ps_invTT_transpose,CC.invTT',prA_TT_Ps_transpose)
            # offers 
            for Ci = 1:J
                prj!B = ws.temps_a[1]
                prj!B .= 0
                for bb=1:nportf
                    if CC.portf[Ci,bb] > 0
                        prj!B[bb] = getenrollmentprob(ui,λ,bb,Ci,CC,ws) #uses temps_j[1], temps_j[2]
                    end
                end
                pr = dot(prA_TT_Ps_invTT_transpose,prj!B)*exp(ws.ll[mm])*wq
                _pr_[Ci] = pr
                _denom_[Ci] += pr
            end
            # fill in outcomes: E(q), E(q+mu), E(v)
            Eq .+= qq .* _pr_
            for j=1:J
                vv = exp(log_σ_e[j])*grid.ep[j,mm]
                rc = view(ws.temps_j[2],1:length(ws.σ_rc))
                rc .= ws.σ_rc .* view(grid.rc,:,mm)
                urc = dot(rc,view(CC.xrc,j,:))
                Ev[j] += (urc+vv)*_pr_[j]
                if !(ttti_admit && JTXpublic[j])
                    zgqj = zγ[j]+s + qq
                    Eshock[j] += (qq + mean(truncated(Normal(),-zgqj,Inf)))*_pr_[j]
                else
                    Eshock[j] += qq*_pr_[j]
                end
            end
        end
    end
    Eq ./= _denom_
    Eshock ./= _denom_
    Ev ./= _denom_
    Eq,Eshock,Ev
end

#for cf's which change weights on acad. indices
function get_γadmit(θstar,θoutcome,w=1.0)
    γadmit0 = θstar[1:J+nz]
    γacademic = θoutcome[2:nz+1]
    γpersonal = γadmit0[J+1:end] - γacademic
    γadmit_new = [γadmit0[1:J];  γpersonal .+ w .* γacademic; ]
    γq = (1-θoutcome[nz+3]) + w*θoutcome[nz+3]
    θout = copy(θstar)
    θout[1:J+nz] = γadmit_new
    θout, γq 
end