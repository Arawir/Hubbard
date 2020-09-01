#ifndef MODEL
#define MODEL

#include "itensor/all.h"
#include "customelectron.h"

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

AutoMPO hubbardHamiltonianAmpo(cElectron &sites,
                       int L, double t, double U, double K, double Jh)
{
    auto ampo = AutoMPO(sites);
    for(int j=1; j<L; j++){
        ampo += +t,"Cdagup",j,"Cup",j+1;
        ampo += +t,"Cdagup",j+1,"Cup",j;
        ampo += +t,"Cdagdn",j,"Cdn",j+1;
        ampo += +t,"Cdagdn",j+1,"Cdn",j;

        ampo += K/2,"S+1",j,"S-1",j+1;
        ampo += K/2,"S-1",j,"S+1",j+1;
        ampo += K,"Sz1",j,"Sz1",j+1;
    }

    if(Args::global().getBool("PBC")){
        ampo += +t,"Cdagup",L,"Cup",1;
        ampo += +t,"Cdagup",1,"Cup",L;
        ampo += +t,"Cdagdn",L,"Cdn",1;
        ampo += +t,"Cdagdn",1,"Cdn",L;

        ampo += K/2,"S+1",L,"S-1",1;
        ampo += K/2,"S-1",L,"S+1",1;
        ampo += K,"Sz1",L,"Sz1",1;
    }

    for(int j=1; j<=L; j++){
        ampo += +U,"Nupdn",j;
        ampo += +Jh,"S01",j;
    }

    return ampo;
}

MPO hubbardHamiltonian(cElectron &sites,
                       int L, double t, double U, double K, double Jh)
{
    return toMPO(hubbardHamiltonianAmpo(sites,L,t,U,K, Jh));
}

MPO generateMz(cElectron &sites, int L)
{
    auto Mz = AutoMPO(sites);

    for(int i=1; i<=L; i++){
        Mz += 0.5*pow(-1,i+1),"Sz",i;
    }

    return toMPO(Mz);
}

double calculateMz(const cElectron &sites, const MPS &psi)
{
    auto Mz = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        Mz += 0.5*pow(-1,i+1),"Sz",i;
    }

    return innerC(psi,toMPO(Mz),psi).real();
}

double calculateDoublon(const cElectron &sites, const MPS &psi)
{
    auto D = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        D += 1.0,"Nupdn",i;
    }

    return innerC(psi,toMPO(D),psi).real();
}
double calculateN(const cElectron &sites, const MPS &psi)
{
    auto D = AutoMPO(sites);

    for(int i=1; i<=psi.length(); i++){
        D += 1.0,"Ntot",i;
    }

    return innerC(psi,toMPO(D),psi).real();
}


#endif // MODEL

