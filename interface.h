#ifndef INTERFACE
#define INTERFACE

#include "itensor/all.h"
#include <ctime>
#include <functional>

double getD(std::string data)
{
    return Args::global().getReal(data);
}

double getR(std::string data)
{
    return Args::global().getReal(data);
}

int getI(std::string data)
{
    return Args::global().getInt(data);
}

bool getB(std::string data)
{
    return Args::global().getBool(data);
}

std::string getS(std::string data)
{
    return Args::global().getString(data);
}

struct Param{
    std::string name;
    std::string type;
    std::string value;
    bool wasSet=false;

    int getInt(){
        if(type!="int"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return atoi(value.c_str());
    }

    double getDouble(){
        if(type!="double"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return atof(value.c_str());
    }

    bool getBool(){
        if(type!="bool"){
            std::cerr << "WARNING: wrong type conversion (" + name + ")!" << std::endl;
        }
        return (bool)atoi(value.c_str());
    }
};

class Parameters{
private:
    std::vector<Param> params;

public:
    Parameters(){ }
    ~Parameters() = default;

    void add(std::string name, std::string type, std::string value){
        params.push_back( Param{name, type, value} );
    }

    void set(int argc, char *argv[]){
        for(int i=1; i<argc; i++){
            set(argv[i]);
        }
        addArgs();
        write();
    }

    bool getBool(std::string name){
        return param(name)->getBool();
    }
    int getInt(std::string name){
        return param(name)->getInt();
    }
    double getDouble(std::string name){
        return param(name)->getDouble();
    }

    void write(){
        std::cout << "---------------------------- System parameters ---------------------------" << std::endl;
        for(auto& param : params){
            std::cout << "  " << param.name << " -> ";
            std::cout << "  " << param.value << " ";
            if(!param.wasSet){ std::cout << "(default)"; }
            std::cout << std::endl;
        }
    }

private:
    void set(std::string data){
        size_t pos = data.find("=");
        param( data.substr(0,pos) )->value = data.substr(pos+1);
        param( data.substr(0,pos) )->wasSet = true;
    }

    void addArgs(){
        for(auto &param : params){
            if(param.type=="bool"){
                Args::global().add(param.name, param.getBool());
            } else if(param.type=="string"){
                Args::global().add(param.name, param.value);
            } else if(param.type=="int"){
                Args::global().add(param.name, param.getInt());
            } else if(param.type=="double"){
                Args::global().add(param.name, param.getDouble());
            }
        }
    }


    Param* param(std::string name){
        for(auto &param : params){
            if(param.name == name) return &param;
        }
        assert("ERROR: Wrong parameter name (" +name + ")!");
        return nullptr;
    }
} Params;



std::vector<std::string> parseInitState(std::string str)
{
    int L = Params.getInt("L");

    std::vector<std::string> state{};
    int i=0;
    state.push_back("");
    while(str[i]!='\0'){
        if(str[i]=='-'){
            state.push_back("");
        } else {
            state.back().append(std::string{str[i]});
        }
        i++;
    }

    if(state.size()==0){
        std::cerr << "  Init state is not delcared!" << std::endl;
    }else if(state.size()<(uint)L){
        uint templateLength = state.size()-1;
        std::cout << "  Init state deduced" << std::endl;
        for(int i=state.size()-1; i<L; i++){
            state.push_back( state.at(i-templateLength) );
        }
    }
    return state;
}

itensor::MPS prepareInitState(Electron &sites)
{
    std::vector<std::string> initState = parseInitState(Args::global().getString("state"));
    int L = Params.getInt("L");

    auto state = itensor::InitState(sites);
    for(int i=1; i<=L; i++){
       state.set(i,initState[i-1]);
    }

    return itensor::MPS(state);
}

struct Experiment
{
    std::string name;
    std::function<void(void)> experiment;

    Experiment(std::string nName) :
        name{nName}
    {

    }

    void operator = (std::function<void(void)> exp)
    {
        experiment = exp;
    }
};

class ExperimentsClass
{
private:
    std::vector<Experiment> experiments;
public:
    ExperimentsClass() :
        experiments {  }
    {

    }
    ~ExperimentsClass() = default;

    Experiment& operator () (std::string name)
    {
        for(auto& exp : experiments){
            if(exp.name == name){ return exp; }
        }
        experiments.push_back( Experiment{name} );
        return experiments.back();
    }

    void run()
    {
        experiment(getS("exp"))->experiment();
    }

    Experiment* experiment(std::string name)
    {
        for(auto& exp : experiments){
            if(exp.name == name){ return &exp; }
        }
        std::cerr << "ERROR: unknown experiment name!" << std::endl;
        assert(false);
        return nullptr;
    }

} Experiments;

class Controller
{
public:
    clock_t  time0, timeLast;
public:

    Controller() :
        time0{ clock() }
      , timeLast{ clock() }
    {

    }
    ~Controller() = default;

    void addPoint(std::string text)
    {
        writePointInfo(text);
    }

private:
    void writePointInfo(std::string text)
    {
        std::cout << "--------------------";
        std::cout << text;
        for(int i=0; i<60-(int)text.size(); i++){
           std::cout << "-";
        }

        std::cout << "(" << (clock()-time0)/(double)CLOCKS_PER_SEC << " [s],"
                  << (clock()-timeLast)/(double)CLOCKS_PER_SEC << " [s])";
        std::cout << std::endl;
        timeLast = clock();
    }
} ExpCon;

#endif // INTERFACE

