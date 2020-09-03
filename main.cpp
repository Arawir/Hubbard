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
            auto sites = Electron(getI("L"));
            auto ampo = AutoMPO(sites);
            double t=getD("t");
            double U=getD("U");
            double J=getD("J");
            double L = getI("L");

            for(int j=1; j<L; j++){
                ampo += -t,"Cdagup",j,"Cup",j+1;
                ampo += -t,"Cdagup",j+1,"Cup",j;
                ampo += -t,"Cdagdn",j,"Cdn",j+1;
                ampo += -t,"Cdagdn",j+1,"Cdn",j;
            }

            for(int j=1; j<=L; j++){
                ampo += +U,"Nupdn",j;
            }

            for(int i = 1; i <= getI("L")-1; ++i){
               ampo += J/2,"S+",i,"S-",i+1;
               ampo += J/2,"S-",i,"S+",i+1;
               ampo +=     J,"Sz",i,"Sz",i+1;

            }
            auto H = toMPO(ampo);
            auto H2 = H;
            printfln("ddddMaximu bond dimension of H is %d",maxLinkDim(H));


            auto psi1 = prepareInitState(sites);
            auto sweeps = prepareSweepClass();


            // start TDVP, either one site or two site algorithm can be used by adjusting the "NumCenter" argument
            println("----------------------------------------GSE-TDVP---------------------------------------");

            std::cout << "  Energy: " << innerC(psi1,H,psi1).real() << std::endl;

            std::cout << " N  " << calculateN(sites,psi1) << " "
                    << " D  " << calculateDoublon(sites,psi1) << " "
                      << " Mz  " << calculateMz(sites,psi1) << std::endl;

            for(double time=0.0; time<=getD("maxtime")+0.001; time+=getD("dtime")){
                if(time<2*getD("dtime")){
                   std::vector<Real> epsilonK = {1E-8, 1E-8,1E-8};
                   addBasis(psi1,H,epsilonK,{"Cutoff",1E-8,"Method","DensityMatrix","KrylovOrd",4,"DoNormalize",true,"Quiet",false});
                }
                std::cout << "  Time: " << time
                          << "  Energy: "
                          << tdvp(psi1,H,im*getD("dtime"),sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",1})
                          << " N " << calculateN(sites,psi1) << " "
                          << " D " << calculateDoublon(sites,psi1) << " "
                          << " Mz  " << calculateMz(sites,psi1) << std::endl;
            }
            printfln("Using overlap = %.10f", real(innerC(psi1,H,psi1)) );
            printfln("Using overlap = %.10f", real(innerC(psi1,H2,psi1)) );
        };



    Params.add("t","double","0.0");
    Params.add("U","double","0.0");
    Params.add("J","double","0.0");

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

