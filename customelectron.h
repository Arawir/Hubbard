#ifndef CUSTOMELECTRON
#define CUSTOMELECTRON

#pragma once

#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

namespace itensor {

class cElectronSite;

using cElectron = BasicSiteSet<cElectronSite>;

class cElectronSite
    {
    Index s;
    public:

    cElectronSite(Index I) : s(I) { }

    cElectronSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,Elec");
        if(args.defined("SiteNumber"))
            {
            ts.addTags("n="+str(args.getInt("SiteNumber")));
            }
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveNf = args.getBool("ConserveNf",conserveQNs);
        auto conserveSz = args.getBool("ConserveSz",conserveQNs);
        if(conserveQNs || conserveNf || conserveSz)
            {
            if(conserveNf && conserveSz)
                {
                s = Index(QN({"Sz", 1},{"Nf",0,-1}),1,
                          QN({"Sz", 2},{"Nf",1,-1}),1,
                          QN({"Sz", 0},{"Nf",1,-1}),1,
                          QN({"Sz", 1},{"Nf",2,-1}),1,
                          QN({"Sz",-1},{"Nf",0,-1}),1,
                          QN({"Sz", 0},{"Nf",1,-1}),1,
                          QN({"Sz",-2},{"Nf",1,-1}),1,
                          QN({"Sz",-1},{"Nf",2,-1}),1,Out,ts);
                }
            else if(conserveNf) // don't conserve Sz
                {
                s = Index(QN({"Nf",0,-1}),1, //Up
                          QN({"Nf",1,-1}),1,
                          QN({"Nf",1,-1}),1,
                          QN({"Nf",2,-1}),1,
                          QN({"Nf",0,-1}),1, //Dn
                          QN({"Nf",1,-1}),1,
                          QN({"Nf",1,-1}),1,
                          QN({"Nf",2,-1}),1,Out,ts);
                }
            else if(conserveSz) //don't conserve Nf, only fermion parity
                {
                s = Index(QN({"Sz", 1},{"Pf",0,-2}),1,
                          QN({"Sz", 2},{"Pf",1,-2}),1,
                          QN({"Sz", 0},{"Pf",1,-2}),1,
                          QN({"Sz", 1},{"Pf",0,-2}),1,
                          QN({"Sz",-1},{"Pf",0,-2}),1,
                          QN({"Sz", 0},{"Pf",1,-2}),1,
                          QN({"Sz",-2},{"Pf",1,-2}),1,
                          QN({"Sz",-1},{"Pf",0,-2}),1,Out,ts);
                }
            else
                {
                s = Index(QN({"Pf",0,-2}),1,
                          QN({"Pf",1,-2}),1,
                          QN({"Pf",1,-2}),1,
                          QN({"Pf",0,-2}),1,
                          QN({"Pf",0,-2}),1,
                          QN({"Pf",1,-2}),1,
                          QN({"Pf",1,-2}),1,
                          QN({"Pf",0,-2}),1,Out,ts);
                }
            }
        else
            {
            s = Index(8,ts);
            }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
    {
        if(state == "0U" || state == "EmpU"){ return s(1); }
        else if(state == "+U" || state == "UpU"){ return s(2);}
        else if(state == "-U" || state == "DnU"){ return s(3);}
        else if(state == "SU" || state == "UpDnU"){ return s(4); }
        else if(state == "0D" || state == "EmpD"){ return s(5);}
        else if(state == "+D" || state == "UpD"){ return s(6);}
        else if(state == "-D" || state == "DnD"){ return s(7);}
        else if(state == "SD" || state == "UpDnD"){ return s(8);}
        else{ Error("State " + state + " not recognized"); }

        return IndexVal{};
    }

        ITensor
        op(std::string const& opname,
           Args const& args) const
        {
        auto sP = prime(s);

        IndexVal EmU(s(1)),
                 EmPU(sP(1)),
                 UpU(s(2)),
                 UpPU(sP(2)),
                 DnU(s(3)),
                 DnPU(sP(3)),
                 UDU(s(4)),
                 UDPU(sP(4)),
                EmD(s(5)),
                EmPD(sP(5)),
                UpD(s(6)),
                UpPD(sP(6)),
                DnD(s(7)),
                DnPD(sP(7)),
                UDD(s(8)),
                UDPD(sP(8));

        ITensor Op(dag(s),sP);

        if(opname == "Nup")
            {
            Op.set(UpU,UpPU,1);
            Op.set(UDU,UDPU,1);
            Op.set(UpD,UpPD,1);
            Op.set(UDD,UDPD,1);
            }
        else
        if(opname == "Ndn")
            {
            Op.set(DnU,DnPU,1);
            Op.set(UDU,UDPU,1);
            Op.set(DnD,DnPD,1);
            Op.set(UDD,UDPD,1);
            }
        else
        if(opname == "Nupdn")
            {
            Op.set(UDU,UDPU,1);
            Op.set(UDD,UDPD,1);
            }
        else
        if(opname == "Ntot")
            {
            Op.set(UpU,UpPU,1);
            Op.set(DnU,DnPU,1);
            Op.set(UDU,UDPU,2);

            Op.set(UpD,UpPD,1);
            Op.set(DnD,DnPD,1);
            Op.set(UDD,UDPD,2);
            }
        else
        if(opname == "Cup")
            {
            Op.set(UpU,EmPU,1);
            Op.set(UDU,DnPU,1);

            Op.set(UpD,EmPD,1);
            Op.set(UDD,DnPD,1);
            }
        else
        if(opname == "Cdagup")
            {
            Op.set(EmU,UpPU,1);
            Op.set(DnU,UDPU,1);

            Op.set(EmD,UpPD,1);
            Op.set(DnD,UDPD,1);
            }
        else
        if(opname == "Cdn")
            {
            Op.set(DnU,EmPU,1);
            Op.set(UDU,UpPU,-1);

            Op.set(DnD,EmPD,1);
            Op.set(UDD,UpPD,-1);
            }
        else
        if(opname == "Cdagdn")
            {
            Op.set(EmU,DnPU,1);
            Op.set(UpU,UDPU,-1);

            Op.set(EmD,DnPD,1);
            Op.set(UpD,UDPD,-1);
            }
        else
        if(opname == "Aup")
            {
            Op.set(UpU,EmPU,1);
            Op.set(UDU,DnPU,1);

            Op.set(UpD,EmPD,1);
            Op.set(UDD,DnPD,1);
            }
        else
        if(opname == "Adagup")
            {
            Op.set(EmU,UpPU,1);
            Op.set(DnU,UDPU,1);

            Op.set(EmD,UpPD,1);
            Op.set(DnD,UDPD,1);
            }
        else
        if(opname == "Adn")
            {
            Op.set(DnU,EmPU,1);
            Op.set(UDU,UpPU,1);

            Op.set(DnD,EmPD,1);
            Op.set(UDD,UpPD,1);
            }
        else
        if(opname == "Adagdn")
            {
            Op.set(EmU,DnPU,1);
            Op.set(UpU,UDPU,1);

            Op.set(EmD,DnPD,1);
            Op.set(UpD,UDPD,1);
            }
        else
        if(opname == "FermiPhase" || opname == "F")
            {
            Op.set(EmU,EmPU,+1);
            Op.set(UpU,UpPU,-1);
            Op.set(DnU,DnPU,-1);
            Op.set(UDU,UDPU,+1);

            Op.set(EmD,EmPD,+1);
            Op.set(UpD,UpPD,-1);
            Op.set(DnD,DnPD,-1);
            Op.set(UDD,UDPD,+1);
            }
        else
        if(opname == "Fup")
            {
            Op.set(EmU,EmPU,+1);
            Op.set(UpU,UpPU,-1);
            Op.set(DnU,DnPU,+1);
            Op.set(UDU,UDPU,-1);

            Op.set(EmD,EmPD,+1);
            Op.set(UpD,UpPD,-1);
            Op.set(DnD,DnPD,+1);
            Op.set(UDD,UDPD,-1);
            }
        else
        if(opname == "Fdn")
            {
            Op.set(EmU,EmPU,+1);
            Op.set(UpU,UpPU,+1);
            Op.set(DnU,DnPU,-1);
            Op.set(UDU,UDPU,-1);

            Op.set(EmD,EmPD,+1);
            Op.set(UpD,UpPD,+1);
            Op.set(DnD,DnPD,-1);
            Op.set(UDD,UDPD,-1);
            }
        else
        if(opname == "Sz0")
            {
            Op.set(UpU,UpPU,+0.5);
            Op.set(DnU,DnPU,-0.5);

            Op.set(UpD,UpPD,+0.5);
            Op.set(DnD,DnPD,-0.5);
            }
        else
        if(opname == "Sz1")
        {
            Op.set(EmU,EmPU,+0.5);
            Op.set(UpU,UpPU,+0.5);
            Op.set(DnU,DnPU,+0.5);
            Op.set(UDU,UDPU,+0.5);

            Op.set(EmD,EmPD,-0.5);
            Op.set(UpD,UpPD,-0.5);
            Op.set(DnD,DnPD,-0.5);
            Op.set(UDD,UDPD,-0.5);
        }
        else
        if(opname == "S01")
        {
            Op.set(UpU,UpPU,+0.25);
            Op.set(DnU,DnPU,-0.25);
            Op.set(UpD,UpPD,-0.25);
            Op.set(DnD,DnPD,+0.25);

            Op.set(UpD,DnPU,0.5);
            Op.set(DnU,UpPD,0.5);
        }
        else
        if(opname == "S+1")
        {
            Op.set(EmD,EmPU,1);
            Op.set(UpD,UpPU,1);
            Op.set(DnD,DnPU,1);
            Op.set(UDD,UDPU,1);
        }
        else
        if(opname == "S-1")
        {
            Op.set(EmU,EmPD,1);
            Op.set(UpU,UpPD,1);
            Op.set(DnU,DnPD,1);
            Op.set(UDU,UDPD,1);
        }
        else
        if(opname == "S+0")
            {
            Op.set(DnU,UpPU,1);

            Op.set(DnD,UpPD,1);
            }
        else
        if(opname == "S-0")
            {
            Op.set(UpU,DnPU,1);

            Op.set(UpD,DnPD,1);
            }
        else
        if(opname == "S20")
            {
            //S dot S on-site
            Op.set(UpU,UpPU,0.75);
            Op.set(DnU,DnPU,0.75);

            Op.set(UpD,UpPD,0.75);
            Op.set(DnD,DnPD,0.75);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

    //
    // Deprecated, for backwards compatibility
    //

    cElectronSite(int n, Args const& args = Args::global())
        {
        *this = cElectronSite({args,"SiteNumber=",n});
        }

    };

//
// Deprecated, for backwards compatability
//

using cHubbardSite = cElectronSite;

using cHubbard = BasicSiteSet<cHubbardSite>;


} //namespace itensor

#endif // CUSTOMELECTRON

