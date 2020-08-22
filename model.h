#ifndef MODEL
#define MODEL

#include "itensor/all.h"

using namespace itensor;


Sweeps prepareSweepClass()
{
    auto sweeps = Sweeps(Args::global().getInt("sweeps"));
    sweeps.maxdim() = Args::global().getInt("maxDim");
    sweeps.mindim() = Args::global().getInt("minDim");
    sweeps.cutoff() = Args::global().getReal("cutoff");
    sweeps.niter() = Args::global().getReal("niter");
    return sweeps;
}

MPO hubbardHamiltonian(Electron &sites,
                       int L, double t, double U)
{
    auto ampo = AutoMPO(sites);
    for(int j=1; j<L; j++){
        ampo += -t,"Cdagup",j,"Cup",j+1;
        ampo += +t,"Cup",j,"Cdagup",j+1;
        ampo += -t,"Cdagdn",j,"Cdn",j+1;
        ampo += +t,"Cdn",j,"Cdagdn",j+1;
    }

    if(Args::global().getBool("PBC")){
        ampo += -t,"Cdagup",L,"Cup",1;
        ampo += +t,"Cup",L,"Cdagup",1;
        ampo += -t,"Cdagdn",L,"Cdn",1;
        ampo += +t,"Cdn",L,"Cdagdn",1;
    }

    for(int j=1; j<=L; j++){
        ampo += +U,"Nupdn",j;
    }

    return toMPO(ampo);
}


double calculateMz(const Electron &sites, const MPS &psi)
{
    auto Mz = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        Mz += 0.5*pow(-1,i+1),"Sz",i;
    }

    return innerC(psi,toMPO(Mz),psi).real();
}

double calculateDoublon(const Electron &sites, const MPS &psi)
{
    auto D = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        D += 1.0,"Nupdn",i;
    }

    return innerC(psi,toMPO(D),psi).real();
}


#endif // MODEL

