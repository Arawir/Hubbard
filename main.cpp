#include "model.h"
#include "interface.h"
#include "tdvp.h"
#include "itensor/all.h"

#include <complex>
#define im std::complex<double>{0.0,1.0}


int main(int argc, char *argv[])
{

    Experiments("DMRG") = [](){
        ExpCon.addPoint("Initialization");

        seedRNG(1);
        auto sites = cElectron( getI("L") );
        auto psi = prepareInitState(sites);
        auto H = hubbardHamiltonian(sites,getI("L"),getD("t"),getD("U"));
        auto sweeps = prepareSweepClass();

        std::cout << "  Energy: " << std::real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  Mz: " << calculateMz(sites,psi) << std::endl;
        std::cout << "  N: " << calculateN(sites,psi) << std::endl;

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        std::cout << "  Energy: " << std::real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  Mz: " << calculateMz(sites,psi) << std::endl;
        std::cout << "  N: " << calculateN(sites,psi) << std::endl;
    };


    Experiments("timeEv") = [](){
        ExpCon.addPoint("Initialization");

        seedRNG(1);
        auto sites = cElectron( getI("L") );
        auto psi = prepareInitState(sites);
        auto H = hubbardHamiltonian(sites,getI("L"),getD("t"),getD("U"));
        auto sweeps = prepareSweepClass();

        ExpCon.addPoint("Time evolution of initial state");

        double time = 0.0;
        while(time <= getD("time")+0.001){
            std::cout << "  StepE ";
            std::cout << time << " ";
            std::cout << calculateMz(sites,psi) << " ";
            std::cout << calculateDoublon(sites,psi) << " ";
            std::cout << std::real(innerC(psi,H,psi)) << std::endl;
            tdvp(psi,H,im*getD("dTime"),sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});

            time += getD("dTime");
        }

        ExpCon.addPoint("Finish");
    };

    Experiments("timeEvBasic") = [](){
        ExpCon.addPoint("Initialization");

        int L = getI("L");
        auto sites = cElectron(L);
        auto psi = prepareInitState(sites);
        auto Hampo = hubbardHamiltonianAmpo(sites,getI("L"),getD("t"),getD("U"));
        auto H = toMPO(Hampo);
        auto Mz = generateMz(sites, L);

        auto expH = toExpH(Hampo,im*getD("dTime"));

        auto args = Args("Method=","DensityMatrix","Cutoff=",1E-9,"MaxDim=",3000);

        for(double time=0.0; time<=getD("time")+0.001; time+=getD("dTime")){
            std::cout << "  StepE ";
            std::cout << time << " ";
            std::cout << std::real(innerC(psi,H,psi)) << " ";
            std::cout << std::real(innerC(psi,Mz,psi)) << " ";
            psi = applyMPO(expH,psi,args);
            psi.noPrime().normalize();
            std::cout << std::endl;
        }

        ExpCon.addPoint("Finish");
    };



    Params.add("t","double","0.0");
    Params.add("U","double","0.0");

    Params.add("L","int","4");
    Params.add("PBC","bool","0");

    Params.add("dTime","double","0.1");
    Params.add("time","double","1.0");

    Params.add("Silent","bool","1");
    Params.add("cutoff","double","1E-6");
    Params.add("sweeps","int","4");
    Params.add("minDim","int","1");
    Params.add("maxDim","int","100");
    Params.add("niter","int","10");
    Params.add("state","string","UpU-DnU");
    Params.add("ConserveNf","bool","0");
    Params.add("ConserveSz","bool","0");
    Params.add("exp","string","1");
    Params.add("ConserveQNs","bool","0");
    Params.add("Method","DensityMatrix","0");

    Params.set(argc,argv);
    Experiments.run();

    return 0;
}

