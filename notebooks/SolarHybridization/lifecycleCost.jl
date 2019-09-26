include("smoothMinMaxAbs.jl")
function lifecycleCost(Q,xL,xU,xts,xa,m)
    #===========================================================================
    # This subroutine calculates the the total lifecycle savings by implementing
    # a solar hybridization strategy: fLCS = NG cost avoided - investment cost
    # Written by: M. D. Stuber, May 17, 2018, rev. May 15, 2019 (Julia 1.0)
    # This code is used in the paper: Stuber (2018), DOI: 10.3390/pr6070076
    # INPUTS
    #   Q: vector of heat outputs from solar array for each time point
    #   X: vector of optimization variable bounds
    #   xts: thermal storage decision variable (h)
    #   xa: aperture area decision variable (m2)
    #   m: capital model to use (=1 for convex, =2 for nonconvex)
    # OUTPUTS
    #   flsc: total lifecycle savings if you implement solar
    ===========================================================================#
    # economic parameters of the model
    r_th = 0.01 # conventional fuel price inflation
    r_d = 0.10 # discount rate
    r_c = 0.065 # capital APR
    term = 10::Int # loan term [years]
    L = 30::Int # total project lifetime [years]
    ng = 7.50/293.1/1.037 # Commercial NG cost $/mmbtu -> $/kWh
    #ng = 3.65/293.1/1.037 # Industrial NG cost $/mmbtu -> $/kWh
    Cts0 = 45.14 #TES cost prefactor
    Ca0  = 425.0 #PTC cost prefactor
    ##################################################
    t = range(1,stop=8760,length=8760)
    qmean = 10.0*1000.0 # mean thermal power demand of IPH [kW]
    qdev = 0.00 # deviation fraction from the mean for ramp up/down (∈[0,1], 0:constant, >0:transient load)
    qp_peak = qmean*(1.0+qdev)# peak IPH demand
    qproc = qmean*(1.0.+qdev*sin.((t.-7.0)/12.0*pi)) # dynamic IPH demand [kW] vector
    qsol = Q*xa/1000.0 # total solar power [kW] vector
    stLoad = qsol-qproc # total storage load [kW] vector
    # initialize arrays
    qstored = zeros(typeof(xts),length(Q)+1) #total thermal energy stored
    qex   = zeros(typeof(xts),length(Q))# excess heat, h(q_s-q_l) in paper
    qloss = copy(qex)# heat loss, q_l in paper
    qanc  = ones(typeof(xts),length(Q))*-qmean# ancillary heat, q_a in paper
    qng = copy(qanc)*-1.0# heat from NG system, q_ng in paper
    dt = 1 # 1h discrete time points, h in paper
    # set up the energy balance equations and calculate SFs
    SFnum=0.0# initialize numerator of solar fraction
    for i in 1:length(Q)
        qex[i]=dt*(stLoad[i])
        qstored[i+1]=minD1(qp_peak*xts,maxD1(qstored[i]+qex[i],0.0))
        qanc[i]=qstored[i]+qex[i]-qstored[i+1]
        qloss[i]=maxD1(0.0,qanc[i])
        qng[i] = -minD1(0.0,qanc[i])
        SFnum += qsol[i]-maxD1(0.0,qstored[i]+qex[i]-qp_peak*xts)
    end
    SF = SFnum/sum(qproc)# smooth solar fraction
    Cp = [ng*(1+r_th)^(i-1) for i=1:L]*sum(qproc) # opex w/o solar (gas cost)
    hL = xL[1] # lower bound on storage (h)
    aL = xL[2] # lower bound on area (m^2)
    hU = xU[1] # upper bound on storage (h)
    aU = xU[2] # upper bound on area (m^2)
    Cts = Cts0*((hL*qp_peak)^(0.91)+((hU*qp_peak)^(0.91)-(hL*qp_peak)^(0.91))/(hU-hL)*(xts-hL))
    Ca = Ca0*((aL)^(0.92)+(aU^(0.92)-aL^(0.92))/(aU-aL)*(xa-aL))
    #
    if m == :convex # convex capital model
        Ccap0 = Cts+Ca # (16) in paper
    else # nonconvex capital model
        Ccap0 = Ca0*xa^(0.92) + Cts0*(xts*qp_peak)^(0.91) # (14) in paper
    end
    # vectorize annual capital costs
    Ccap = [r_c*Ccap0*(1.0+r_c/12.0)^(12*term)/((1+r_c/12)^(12*term)-1) for i=1:L]
    if term<L Ccap[term+1:L,1] = zeros(L-term,1)  end # set capital cost (debt service) =0 after financing is paid
    disc = [(1+r_d)^(-i) for i = 1:L]# vectorize discount factor for time value
    flcs = (SF*Cp - Ccap)'*disc # (13) in paper
    return flcs
end
