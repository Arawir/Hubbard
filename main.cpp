#include "model.h"
#include "interface.h"
#include "tdvp.h"
#include "basisextension.h"

#include <complex>
#define im std::complex<double>{0.0,1.0}


int main(int argc, char *argv[])
{
    Experiments("dmrg") = [](){

        ExpCon.addPoint("Initialization");
        auto [sites, psi, H, sweeps] = prepareExpBasic();
        ExpCon.setSites(sites); ExpCon("E") = H;
        ExpCon.calc(psi, oMode::b, "E","Mz/L","N");

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        ExpCon.calc(psi, oMode::b, "E","Mz/L","N");
    };

    Experiments("timeEv") = [](){
        ExpCon.addPoint("Initialization");

        auto sites = Electron(getI("L"));
        auto H = hubbardHamiltonian(sites,getI("L"),getD("t"),getD("U"));
        auto Hdist = hubbardHamiltonianWithDist(sites,getI("L"),getD("t"),getD("U"));
        auto psi = prepareInitState(sites);
        auto sweeps = prepareSweepClass();

        ExpCon.setSites(sites); ExpCon("E") = H;

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,Hdist,sweeps);

        ExpCon.addPoint("Starting TDVP");

        double energy = innerC(psi,H,psi).real();

        for(double time=0; time<=getD("maxtime")+getD("dtime")+0.001; time+=getD("dtime")){
            ExpCon.calc(psi,oMode::b,"t:",time,"E","N","D3","maxDim:",maxLinkDim(psi),"Ni:",oMode::a,"N1:L");

            if(time<2*getD("dtime")){
               std::vector<Real> epsilonK = {getD("cutoff"),getD("cutoff"),getD("cutoff")};
               addBasis(psi,H,epsilonK,{"Cutoff",getD("cutoff"),"Method","DensityMatrix","KrylovOrd",4,"DoNormalize",true,"Quiet",getB("Silent")});
            }
            energy = tdvp(psi,H,im*getD("dtime"),sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});
        }
        ExpCon.addPoint("Finish");
    };



    Params.add("t","double","0.0");
    Params.add("U","double","0.0");

    Params.add("L","int","4");
    Params.add("PBC","bool","0");

    Params.add("dtime","double","0.1");
    Params.add("maxtime","double","1.0");

    Params.add("Silent","bool","1");
    Params.add("cutoff","double","1E-6");
    Params.add("sweeps","int","4");
    Params.add("minDim","int","1");
    Params.add("maxDim","int","100");
    Params.add("niter","int","10");
    Params.add("state","string","Up-Dn");
    Params.add("ConserveNf","bool","0");
    Params.add("ConserveSz","bool","0");
    Params.add("exp","string","1");
    Params.add("ConserveQNs","bool","0");

    Params.set(argc,argv);
    prepareObservables();
    Experiments.run();

    return 0;
}
