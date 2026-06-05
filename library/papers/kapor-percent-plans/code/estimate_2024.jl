
function obj_new(θ,CC=CC,popwt_wv2=popwt_wv2,IW=1:II,likwt=ones(II),aid_data = Xj.avgfinaid_j,W=100 .* ones(J))
    println("θ=$(round.(ForwardDiff.value.(θ),digits=6))")
    II = length(CC.ttt)
    lik = zeros(eltype(θ),II)
    aid = zeros(eltype(θ),II)
    Threads.@threads for ii=1:II
        (lik[ii], aid[ii]) = obj_i_new(θ,CC,ii,IW[ii])
    end
    obj = -dot(lik,likwt)
    for j=1:J
        inds_j = (CC.Ci .== j)
        aid_model = dot(aid[inds_j], popwt_wv2[inds_j])/sum(popwt_wv2[inds_j])
        moment = W[j]*(aid_model - aid_data[j])^2
        #println("j=$j: aid_model=$(ForwardDiff.value(aid_model)), aid_data=$(aid_data[j]), g=$(ForwardDiff.value(moment))")
        obj += moment
    end
    println("obj=$(ForwardDiff.value(obj))")
    obj
end

function obj_i_new(θ,CC,ii,iw=ii)
    T = eltype(θ)
    ws::WS{T} = get!(WSS[iw],T) do 
        WS(CC,iw,T)
    end
    #unpack parameters
    γadmit = view(θ,1:J+nz)  #first j elements are cutoffs; next nz are zi 
    γs = SVector{3}(θ[J+nz+1],θ[J+nz+2],θ[J+nz+3])
    γq!s = SVector{3}(θ[J+nz+4],θ[J+nz+5],θ[J+nz+6])
    λ = 1/(1+exp(-θ[J+nz+7])) #correlation parameter for final-period shocks
    γcfixed = SVector{3}(θ[J+nz+8],θ[J+nz+9],θ[J+nz+10]) #app costs
    γcvariable = SVector{3}(θ[J+nz+11],θ[J+nz+12],θ[J+nz+13])
    #σAppShocks = exp(θ[J+nz+14]) #"smoothing" parameter for application costs
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
    sig_q!s = sqrt(var_q!s)
    c0 = dot(γcfixed, view(CC.zi_info,ii,:))
    c1 = dot(γcvariable, view(CC.zi_info,ii,:))
    zγ = ws.zγ; mul!(zγ,CC.zi[ii],γadmit)  #all schools
    #compute probability of observed B | observed A, q (ignoring binomial terms, which don't depend on parameters)
    zγA = ws.bvec
    inds = CC.inv_inds_iA[ii]
    b_obs = view(CC.biA,inds)
    mul!(zγA,view(CC.ziA,inds,:),γadmit)
    #loop over draws mm=1:MM (s,y,ε^a,awareness), compute likelihood
    grid = CC.grid[ii]
    for mm=1:ndraws
        s = grid.s[mm]*sqrt(var_s)
        ui, aid_m = getu(βx,ws.σ_rc,log_σ_e,βp,coef_aidamount,CC,ii,mm,ws)
        ws.UAd[mm] .= getvaladmit(ui,λ,ws,CC)
        mul!(ws.invT_UAd[mm],CC.invTT,ws.UAd[mm]) #expected payoff of offer, needed for app likelihood
        #error("this fuckin line above")
        #likelihood of awareness draws
        lik_aware = getFinaidAwarenessLikelihood(view(grid.aware,:,mm),βx_aware,βy_aware,σ_βi0_aware,CC.xij[ii],grid.y[mm],grid.efc[mm],ws.temps_j,nothing)
        ws.ll[mm] = log(lik_aware) - log(grid.pr0_aware[mm])
        #app decision
        Pq = get_Pq(s,sig_q!s,zγ,ws,ii,CC) #ws.Pq
        V = getvalapp(ws.invT_UAd[mm],Pq,ws,CC,mm)
        V .-= getappcosts(c0,c1,CC,ii,ws) #deterministic part of costs
        ws.ValApp_allM[mm] .= V
        ll_a = logliklogit(V,CC.Aij[ii],σAppShocks) #likelihood of application, draw mm
        ws.ll[mm] += ll_a
        #likelihood of admission
        ws.logprB!q .= 0
        lik_b = zero(eltype(θ))
        for (ind_q,qq) in enumerate(Qnodes)
            for j=1:length(inds)
                (CC.ttt[ii] && JTXpublic[j]) && continue
                zgqj = zγA[j]+s + qq*sqrt(2)*sig_q!s
                ws.logprB!q[ind_q] += (b_obs[j] ? logccdf(Normal(),-zgqj) : logcdf(Normal(),-zgqj))
            end   
            lik_b += exp(ws.logprB!q[ind_q])*Qweights[ind_q]    
        end
        ws.ll[mm] += log(lik_b)
        #final matriculation decision
        if CC.wv2[ii]
            ws.ll[mm] += log( getenrollmentprob(ui,λ,CC.Bij[ii],CC.Ci[ii],CC,ws))
            CC.Ci[ii] > 0 && (ws.aid[mm] = aid_m[CC.Ci[ii]])
        end
    end
    loglik,avgaid = pack_ll_aid(ws,CC,ii)
end

@inline function pack_ll_aid(ws,CC,ii)
    T = eltype(ws.ll)
    ll = ws.ll
    aid = ws.aid
    lmax = maximum(ll)
    expll = ws.temps_m[1]
    sum_expll = zero(T)
    for mm=1:ndraws
        expll[mm] = exp(ll[mm]-lmax)
        sum_expll += expll[mm]
    end
    avgaid = zero(T)
    if CC.wv2[ii] && (CC.Ci[ii] > 0)
        for mm=1:ndraws
            avgaid += aid[mm]*expll[mm] / sum_expll
        end
    end
    loglik = log(sum_expll) + lmax
    return loglik,avgaid
end

@inline function get_Pq(s,sig_q!s,zγ,ws,ii,CC,TTTi = CC.ttt[ii], wq=1.0)
    T = promote_type(typeof(s),typeof(sig_q!s),eltype(zγ))
    for (ind_q,qq) in enumerate(Qnodes)
        #fill in Pq, prob of admission to all
        Pq_n = ws.Pq[ind_q]
        Pq_n .= 1
        @inbounds for j=1:J
            #(JTXpublic[j] && CC.ttt[ii]) && continue
            (JTXpublic[j] && TTTi) && continue
            zgqj = zγ[j]+s + qq*sqrt(2)*sig_q!s*wq
            prj = ccdf(Normal(),-zgqj)
            for aa = 1:size(Pq_n,1)
                Pqnaj = one(T)
                nj = CC.portf_T[aa,j]
                for counter = nj : -1 : 1
                    Pqnaj *= prj
                end
                Pq_n[aa] *= Pqnaj
            end
        end
    end
    ws.Pq
end

function getu(βx,σ_rc,log_σ_e,βp,coef_aidamount,CC,ii,mm,ws)
    grid = CC.grid[ii]
    u = ws.u
    aid_m = ws.aid_m
    aid_m .= max.(zero(eltype(u)), (CC.listprice .- grid.efc[mm]).* exp.(coef_aidamount))
    aid_m .*= view(grid.aware,:,mm)
    netprice = ws.temps_j[1]
    netprice .= CC.listprice .- aid_m
    pricecoef = -log(1+exp(βp[1] + βp[2]*grid.y[mm]))
    rc = view(ws.temps_j[2],1:length(σ_rc))
    rc .= σ_rc .* view(grid.rc,:,mm)
    σ_e = ws.temps_j[3]
    σ_e .= exp.(log_σ_e)
    mul!(u,CC.xij[ii],βx)
    mul!(u,CC.xrc,rc,1,1)
    u .+= netprice .* pricecoef
    u .+= σ_e .* view(grid.ep,:,mm)
    u,aid_m
end

function getenrollmentprob(ui,λ,Bij::Int,Ci::Int,CC,ws)
    T = promote_type(eltype(ui),typeof(λ))
    if Bij==1
        return one(T)
    else
        umax = maximum(ui)
        expλu = ws.temps_j[1]
        expλu .= exp.((ui .- umax)./λ)
        inner_terms = ws.temps_j[2]
        for jj=1:J
            inner_terms[jj] = CC.portf[jj,Bij] * expλu[jj]
        end
        sum_inner_terms = sum(inner_terms)
        _G = sum_inner_terms^λ + exp(-umax)
        if Ci > 0
            return inner_terms[Ci] * sum_inner_terms^(λ - 1) / _G
        else
            return exp(-umax) / _G
        end
    end
end

function logliklogit(v,Ai,σ)
    maxv = maximum(v)
    denom = zero(eltype(v))
    @inbounds @simd for t=1:length(v)
        denom += exp((v[t]-maxv)/σ)
    end
    s = (v[Ai]-maxv)/σ - log(denom)
end

function getappcosts(c0, c1, CC, ii, ws)
    tmp = ws.temps_j[1]
    tmp .= c1
    out = ws.temps_a[1]
    mul!(out,CC.portf',tmp)
    for aa=2:length(out)
        out[aa] += c0
    end
    out
end

function getvalapp(invT_UAd,Pq,ws,CC,mm)
    #inplace V = CC.TT*Diagonal(P)*CC.invTT*UAd
    Ps = ws.Ps[mm]
    Ps .= 0
    for ind_q = 1:length(Qnodes)
        Ps .+= Pq[ind_q].*Qweights[ind_q]
    end
    tmp = ws.temps_a[2]
    tmp .= Ps .* invT_UAd
    V = ws.ValApp
    mul!(V,CC.TT,tmp)
    V
end

# function getpradmit!s(s,sig_q,sig_q!s,Pq,ws)
#     #∫P(q) f(q|s) dq ≈ Σ_n w_n P(q_n)f(q_n|s)/f(q_n)
#     Ps = ws.temps_a[1]
#     Ps .= 0
#     denom = 0.0
#     for ind_q = 1:length(Qnodes)
#         qq = Qnodes[ind_q]*sqrt(2)*sig_q
#         ff = pdf(Normal(s,sig_q!s),qq)/pdf(Normal(0.0,sig_q),qq)
#         ww = Qweights[ind_q]
#         denom += ff * ww
#         pq = Pq[ind_q]
#         Ps .+= pq .* ww .* ff
#     end
#     Ps ./= denom
# end

function getprq!s(prq!s, s,sig_q,sig_q!s)
    denom = zero(typeof(s))
    for ind_q = 1:length(Qnodes)
        qq = Qnodes[ind_q]*sqrt(2)*sig_q
        ff = pdf(Normal(s,sig_q!s),qq)/pdf(Normal(0.0,sig_q),qq)
        ww = Qweights[ind_q]
        denom += ff * ww
        prq!s[ind_q] = ww*ff
    end
    prq!s ./= denom
end
    
#value of offer sets B before learning final shocks
#fast version; returns ws.avec1
function getvaladmit(ui,λ,ws,CC)
    umax = maximum(ui)
    expλu = ws.temps_j[2]
    expλu .= exp.((ui .- umax)./λ )
    G = ws.temps_a[1]
    mul!(G,CC.portf',expλu)
    G .= exp.(λ .* log.(G)) #G.^λ
    G .+= exp(-umax)
    UAd = ws.temps_a[2] #view(ws.UAd,:,mm)
    UAd .= umax .+ log.(G)
    UAd[1] = 0
    UAd
end

function getFinaidAwarenessLikelihood(aware,βx_aware,βy_aware,σ_βi0_aware,xij,ym,efcm,temps_j::Array{Array{T,1},1}, portf_A, skip_known=false, skip_unknown=false,listprice=listprice) where T
    tmp1::Array{T,1} = temps_j[1]
    tmp2::Array{T,1} = temps_j[2]
    lik = 0.0
    x_terms::Array{T,1} = temps_j[3]; mul!(x_terms,xij,βx_aware)
    for (qq,ww) in zip(Qnodes,Qweights)
        βi0aware = qq*sqrt(2)*σ_βi0_aware
        tmp1 .= x_terms .+ βy_aware*ym .+ βi0aware
        tmp2 .= logistic.(tmp1)
        ll_qq = 0.0
        for jj=1:J
            skip_known && (portf_A[jj]>0) && continue
            skip_unknown && (portf_A[jj]==0) && continue
            (efcm > listprice[jj]) && continue #no likelihood contribution if not eligible for aid
            if aware[jj]
                ll_qq += log(tmp2[jj])
            else
                ll_qq += log(1-tmp2[jj])
            end
        end
        lik += exp(ll_qq) * ww
    end
    lik
end

# @inline function getFinaidAwarenessProbs(βx_aware,βy_aware,σ_βi0_aware,xij,ym,temps_j::Array{Array{T,1},1}) where T
#     tmp1 = temps_j[1]
#     tmp2 = temps_j[2]
#     x_terms = temps_j[3]; mul!(x_terms,xij,βx_aware)
#     pr = temps_j[4]; pr .= 0 #likelihood to fill in
#     for (qq,ww) in zip(Qnodes,Qweights)
#         βi0aware = qq*sqrt(2)*σ_βi0_aware
#         tmp1 .= x_terms .+ βy_aware*ym .+ βi0aware
#         tmp2 .= logistic.(tmp1)
#         pr .+= tmp2 .* ww
#     end
#     pr
# end

function getFinaidAwarenessDraws!(βx_aware,βy_aware,σ_βi0_aware,CC,ii::Int,ws=WSS[ii][Float64])
    Aij = view(CC.portf,:,CC.Aij[ii])
    FinaidAij = view(CC.portf,:,CC.FinaidAij[ii])
    grid = CC.grid[ii]
    grid.pr0_aware .= 1
    for mm=1:ndraws
        aid_eligible = CC.listprice .- CC.grid[ii].efc[mm] .> 0
        βi0_aware = randn()*σ_βi0_aware
        x_terms = CC.xij[ii]*βx_aware
        pr0 = logistic.( x_terms .+ βy_aware*grid.y[mm] .+ βi0_aware)
        for jj=1:J
            if Aij[jj]>0
                aware_jm = FinaidAij[jj]>0 #observe awareness in this case
            else
                aware_jm = rand(Bernoulli(pr0[jj]))
            end
            grid.aware[jj,mm] = aware_jm
        end
        grid.pr0_aware[mm] = getFinaidAwarenessLikelihood(grid.aware[:,mm],βx_aware,βy_aware,σ_βi0_aware,CC.xij[ii],grid.y[mm],grid.efc[mm],ws.temps_j, Aij, true)
    end
end

#convenience function!
function getFinaidAwarenessDraws!(θ,CC,IW=1:II)
    βx_aware = view(θ,J+nz+17+nx+nrc:J+nz+16+2nx+nrc)
    βy_aware = θ[J+nz+17+2nx+nrc]
    σ_βi0_aware = exp(θ[J+nz+18+2nx+nrc])
    for ii=1:length(CC.ttt)
        iw = IW[ii]
        getFinaidAwarenessDraws!(βx_aware,βy_aware,σ_βi0_aware,CC,ii,WSS[iw][Float64])
    end
end

