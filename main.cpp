#include "model.h"
#include "interface.h"
#include "tdvp.h"
#include "basisextension.h"


void tdvpStepWithBasisExtensionIfNeeded(MPS &psi, MPO &H, double dTime, Sweeps &sweeps)
{
    if(maxLinkDim(psi)<10){
       std::vector<Real> epsilonK = {getD("cutoff"),getD("cutoff"),getD("cutoff")};
       addBasis(psi,H,epsilonK,{"Cutoff",getD("cutoff"),"Method","DensityMatrix","KrylovOrd",4,"DoNormalize",true,"Quiet",true});
    }
    tdvp(psi,H,im*dTime,sweeps,{"DoNormalize",true,"Quiet",true,"NumCenter",2});
}

int main(int argc, char *argv[])
{
    Experiments("dmrg") = [](){
        auto [sites, psi, H, sweeps] = prepareExpBasic();
        ExpCon.setSites(sites); ExpCon("E") = H;
        ExpCon.calc(psi, oMode::b, "E","Mz/L","N");

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,H,sweeps);

        ExpCon.addPoint("Output data");
        ExpCon.calc(psi, oMode::b, "E","Mz/L","N");
    };

    Experiments("timeEv") = [](){
        auto [sites, psi, H, sweeps] = prepareExpBasic();
        auto Hdist = hubbardHamiltonianWithDist(sites,getI("L"),getD("t"),getD("U"));
        ExpCon.setSites(sites); ExpCon("E") = H;

        ExpCon.addPoint("Starting DMRG");
        dmrg(psi,Hdist,sweeps);

        ExpCon.addPoint("Starting TDVP");
        for(double time=0; time<=getD("maxtime")+getD("dtime")+0.001; time+=getD("dtime")){
            ExpCon.calc(psi,oMode::b,"t:",time,"rtime","mem","E","N","D3","dim","Ni:",oMode::a,"N1:L");
            tdvpStepWithBasisExtensionIfNeeded(psi,H,getD("dtime"),sweeps);
        }
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

    Params.add("PBSenable","bool","0");
    Params.add("PBSjobid","int","0");

    Params.set(argc,argv);
    prepareObservables();
    Experiments.run();

    return 0;
}
