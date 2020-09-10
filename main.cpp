#include "model.h"
#include "interface.h"
#include "tdvp.h"
#include "basisextension.h"

#include <complex>
#define im std::complex<double>{0.0,1.0}


int main(int argc, char *argv[])
{

    Experiments("DMRG") = [](){
        ExpCon.addPoint("Initialization");

        seedRNG(1);
        auto sites = Electron( getI("L") );
        auto psi = prepareInitState(sites);
        auto H = hubbardHamiltonian(sites,getI("L"),getD("t"),getD("U"));
        auto sweeps = prepareSweepClass();

        std::cout << "  Energy: " << std::real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  Mz: " << calculateMzPerL(sites,psi) << std::endl;

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        std::cout << "  Energy: " << std::real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  Mz: " << calculateMzPerL(sites,psi) << std::endl;
    };
    Experiments("timeEv") = [](){
        ExpCon.addPoint("Initialization");

        auto sites = Electron(getI("L"));
        auto H = hubbardHamiltonian(sites,getI("L"),getD("t"),getD("U"));
        auto psi = prepareInitState(sites);
        auto sweeps = prepareSweepClass();

        ExpCon.addPoint("Starting TDVP");

        double energy = innerC(psi,H,psi).real();

        for(double time=0; time<=getD("maxtime")+getD("dtime")+0.001; time+=getD("dtime")){
            std::cout << "  t: " << time << " ";
            std::cout << energy << " ";
            std::cout << innerC(psi,H,psi).real() << " ";
            std::cout << calculateNPerL(sites,psi) << " ";
            std::cout << calculateDPerL(sites,psi) << " ";
            std::cout << calculateMsPerL(sites,psi) << " ";
            std::cout << std::endl;

            if(time<2*getD("dtime")){
               std::vector<Real> epsilonK = {getD("cutoff"),getD("cutoff"),getD("cutoff")};
               addBasis(psi,H,epsilonK,{"Cutoff",getD("cutoff"),"Method","DensityMatrix","KrylovOrd",4,"DoNormalize",true,"Quiet",getB("Silent")});
            }
            energy = tdvp(psi,H,im*getD("dtime"),sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",1});
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
    Experiments.run();

    return 0;
}
