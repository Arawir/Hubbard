#ifndef MODEL
#define MODEL

#include "itensor/all.h"
#include "interface.h"

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
MPO hubbardHamiltonianWithDist(Electron &sites,
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

    for(int j=11; j<=L; j++){
        ampo += t*1E6,"Ntot",j;
    }

    return toMPO(ampo);
}

void prepareObservables()
{
    ExpCon("Mz/L") = [](const Electron &sites){
        auto ampo = AutoMPO(sites);

        for(int i=1; i<=sites.length(); i++){
            ampo += 0.5/sites.length()*pow(-1,i+1),"Sz",i;
        }

        return toMPO(ampo);
    };
    ExpCon("Ms/L") = [](const Electron &sites){
        auto ampo = AutoMPO(sites);

        for(int i=1; i<=sites.length(); i++){
            ampo += +1.0/2.0/sites.length()*pow(-1,i),"Nup",i;
            ampo += -1.0/2.0/sites.length()*pow(-1,i),"Ndn",i;
        }

        return toMPO(ampo);
    };
    ExpCon("D/L") = [](const Electron &sites){
        auto ampo = AutoMPO(sites);

        for(int i=1; i<=sites.length(); i++){
            ampo += 1.0/sites.length(),"Nupdn",i;
        }

        return toMPO(ampo);
    };
    ExpCon("N") = [](const Electron &sites){
        auto ampo = AutoMPO(sites);

        for(int i=1; i<=sites.length(); i++){
            ampo += 1.0,"Ntot",i;
        }

        return toMPO(ampo);
    };
    ExpCon("D3") = [](const Electron &sites){
        auto ampo = AutoMPO(sites);

        ampo += 1.0/3.0,"Nupdn",1;
        ampo += 1.0/3.0,"Nupdn",2;
        ampo += 1.0/3.0,"Nupdn",3;

        return toMPO(ampo);
    };
    ExpCon("N1:L") = [](const Electron &sites){
        std::vector<MPO> out;

        for(int i=1; i<=sites.length(); i++){
            auto ampo = AutoMPO(sites);
            ampo += 1.0,"Ntot",i;
            out.push_back( toMPO(ampo) );
        }

        return out;
    };
}

std::tuple<Electron,MPS,MPO,Sweeps> prepareExpBasic()
{
    seedRNG(1);
    auto sites = Electron( getI("L") );
    auto psi = prepareInitState(sites);
    auto H =hubbardHamiltonian(sites,getI("L"),getD("t"),getD("U"));
    auto sweeps = prepareSweepClass();

    return std::make_tuple( sites,psi,H,sweeps );
}



#endif // MODEL

