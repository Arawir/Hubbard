#include "model.h"
#include "interface.h"
#include "tdvp.h"
#include "itensor/all.h"
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
        std::cout << "  Mz: " << calculateMz(sites,psi) << std::endl;

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        std::cout << "  Energy: " << std::real(innerC(psi,H,psi)) << std::endl;
        std::cout << "  Mz: " << calculateMz(sites,psi) << std::endl;
    };


    Experiments("timeEv") = [](){
        ExpCon.addPoint("Initialization");

        seedRNG(1);
        auto sites = Electron( getI("L") );
        auto psi = prepareInitState(sites);
        auto H = hubbardHamiltonian(sites,getI("L"),getD("t"),getD("U"));
        auto sweeps = prepareSweepClass();

        ExpCon.addPoint("Time evolution of initial state");

        double time = 0.0;
        while(time <= getD("maxtime")+0.001){
            std::cout << "  StepE ";
            std::cout << time << " ";
            std::cout << calculateMz(sites,psi) << " ";
            std::cout << calculateDoublon(sites,psi) << " ";
            std::cout << std::real(innerC(psi,H,psi)) << std::endl;
            tdvp(psi,H,im*getD("dtime"),sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});

            time += getD("dtime");
        }

        ExpCon.addPoint("Finish");
    };

    Experiments("timeEv3") = [](){


            auto sites = SpinHalf(getI("L")*2);
            auto ampo = AutoMPO(sites);

            for(int i = 1; i <= getI("L")*2-2; ++ i){
               ampo += 0.5,"S+",i,"S-",i+2;
               ampo += 0.5,"S-",i,"S+",i+2;
               ampo +=     "Sz",i,"Sz",i+2;
            }
            auto H = toMPO(ampo);
            printfln("Maximum bond dimension of H is %d",maxLinkDim(H));


            auto state = InitState(sites);
            state.set(1,"Up");
            for(int i = 2; i < getI("L")*2; i=i+2){
               if((i/2)%2==1){
                   state.set(i,"Dn");
                   state.set(i+1,"Dn");
                } else {
                   state.set(i,"Up");
                   state.set(i+1,"Up");
                }
            }
            state.set(getI("L")*2,"Up");

            auto psi1 = MPS(state);



            auto sweeps = prepareSweepClass();


            // start TDVP, either one site or two site algorithm can be used by adjusting the "NumCenter" argument
            println("----------------------------------------GSE-TDVP---------------------------------------");

            std::cout << "  Energy: " << innerC(psi1,H,psi1).real() << std::endl;


            for(double time=0.0; time<=getD("maxtime")+0.001; time+=getD("dtime")){
                if(time<2*getD("dtime")){
                   std::vector<Real> epsilonK = {1E-12, 1E-12};
                   addBasis(psi1,H,epsilonK,{"Cutoff",1E-12,"Method","DensityMatrix","KrylovOrd",3,"DoNormalize",true,"Quiet",true});
                }
                std::cout << "  Time: " << time
                          << "  Energy: "
                          << tdvp(psi1,H,im*getD("dtime"),sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",1})
                          << std::endl;
            }
            printfln("Using overlap = %.10f", real(innerC(psi1,H,psi1)) );
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

