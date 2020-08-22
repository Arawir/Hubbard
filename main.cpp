#include "model.h"
#include "interface.h"
#include "tdvp.h"

#include <complex>

#define im std::complex<double>{0.0,1.0}

void exp1()
{
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
}

void exp2()
{
    printfln("------------------------------ Initial State -----------------------------");
    clock_t  timecpu = clock();

    seedRNG(1);
    double t = Args::global().getReal("t");
    double U = Args::global().getReal("U");
    int L = Args::global().getInt("L");

    double dTime = Args::global().getReal("dTime");
    double maxTime = Args::global().getReal("time");

    auto sites = Electron( Args::global().getInt("L") );
    auto psi = prepareInitState(sites);
    auto H = hubbardHamiltonian(sites,L,t,U);
    auto sweeps = prepareSweepClass();

    std::cout << "  Energy: " << std::real(innerC(psi,H,psi)) << std::endl;
    std::cout << "  Mz: " << calculateMz(sites,psi) << std::endl;
    std::cout << "  Nupdn: " << calculateDoublon(sites,psi) << std::endl;
    std::cout << "  Time: " << (clock()-timecpu)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;

    printfln("------------------ Time evolution of initial state ----------------");
    timecpu = clock();

    double time = 0.0;
    while(time <= maxTime+0.001){
        std::cout << time << " ";
        std::cout << calculateMz(sites,psi) << " ";
        std::cout << calculateDoublon(sites,psi) << std::endl;
        tdvp(psi,H,im*dTime,sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",getI("L")/2});

        time += dTime;
    }

    std::cout << "  Energy after TVDP: " << std::real(innerC(psi,H,psi))  << std::endl;
    std::cout << "  Time: " << (clock()-timecpu)/(double)CLOCKS_PER_SEC << " [s]" << std::endl;
}

void run()
{
    switch (getI("exp")){
    case 1 :
        exp1();
        break;
    case 2 :
        exp2();
        break;
    default :
        std::cerr << "ERROR: epxeriment was not selected!" << std::endl;
        break;
    }
}


int main(int argc, char *argv[])
{
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
    Params.add("state","string","Up-Dn");
    Params.add("ConserveNf","bool","0");
    Params.add("ConserveSz","bool","0");
    Params.add("exp","int","1");
    Params.add("ConserveQNs","bool","0");

    Params.set(argc,argv);

    run();

    return 0;
}

