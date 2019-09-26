include("smoothMinMaxAbs.jl")
function iphProcessSmooth(Q,xts,xa)
    #===========================================================================
    # This subroutine calculates the (smooth) solar fraction as a stand-alone
    # for use as a constraint in the optimization problem.  The solar fraction
    # calculation is identical to that in lifecycleCost.jl.
    # Written by: M. D. Stuber, Apr 14, 2018, rev. May 15, 2019
    # This code is used in the paper: Stuber (2018), DOI: 10.3390/pr6070076
    # INPUTS
    #   Q vector of heat outputs from solar array for each time point
    #   xts thermal storage decision variable (h)
    #   xa aperture area decision variable (m2)
    # OUTPUTS
    #   SF solar fraction
    ===========================================================================#
    t = range(1,stop=8760,length=8760)
    qmean = 10.0*1000.0 # mean thermal power demand of IPH [kW]
    qdev = 0.0 # deviation fraction from the mean for ramp up/down
    qp_peak = qmean*(1.0+qdev)
    qproc = qmean*(1.0.+qdev*sin.((t.-7.0)/12.0*pi)) # dynamic IPH demand [kW]
    qsol = Q*xa/1000.0                               # total solar power [kW]
    stLoad = qsol-qproc                              # total storage load [kW]
    qstored = zeros(typeof(xts),length(Q)+1)        #total thermal energy stored
    qex   = zeros(typeof(xts),length(Q))            #total excess thermal energy
    qloss = copy(qex)                               #unutilazable thermal energy
    qanc  = ones(typeof(xts),length(Q))*-qmean      #total ancillary energy needed
    qng = copy(qanc)*-1.0
    dt = 1                                          # 1h discrete time points
    SFnum = 0.0                                     #initialize solar fraction
    # set up the energy balance equations
    for i in 1:length(Q)
        qex[i]=dt*(stLoad[i])
        qstored[i+1]=minD1(qp_peak*xts,maxD1(qstored[i]+qex[i],0.0))
        qanc[i]=qstored[i]+qex[i]-qstored[i+1]
        qloss[i]= maxD1(0,qanc[i])
        qng[i] = -minD1(0,qanc[i])
        SFnum += qsol[i]-maxD1(0.0,qstored[i]+qex[i]-qp_peak*xts)
    end
    SF = SFnum/sum(qproc)
    return SF
end
