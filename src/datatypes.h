#ifndef DATATYPES_H
#define DATATYPES_H

#include <string>
#include <vector>
#include <sstream>
#include <map>
#include "node.h"

struct DataInfo {
public:
    double sc = 1.0;
    
    struct DiffStruct{
        double k0 = 2.0e4, Q = 100.0;
    };

    struct ScaleDefault{
        double val = 1.0;
    };

    struct diPair {
        std::string str;
        int i;
    };

    std::map<std::string, std::string>  nameMap;
    std::vector<diPair>  pressureMap;
    std::map<std::string, double> molMass, add;
    std::map<std::string, ScaleDefault> scale;
    std::map<std::string, DiffStruct> diff;
    std::vector<std::string> GemsStatus;

    DataInfo() {
        GemsStatus.push_back("NO_GEM_SOLVER");
        GemsStatus.push_back("NEED_GEM_AIA");
        GemsStatus.push_back("OK_GEM_AIA");
        GemsStatus.push_back("BAD_GEM_AIA");
        GemsStatus.push_back("ERR_GEM_AIA");
        GemsStatus.push_back("NEED_GEM_SIA");
        GemsStatus.push_back("OK_GEM_SIA");
        GemsStatus.push_back("BAD_GEM_SIA");
        GemsStatus.push_back("ERR_GEM_SIA");
        GemsStatus.push_back("T_ERROR_GEM");

        // from MELCOR to GEMS
        nameMap["XE"]  = "Xe";
        nameMap["CS"]  = "Cs";
        nameMap["BA"]  = "Ba";
        nameMap["I2"]  = "I2";
        nameMap["TE"]  = "Te";
        nameMap["RU"]  = "Ru";
        nameMap["MO"]  = "Mo";
        nameMap["CE"]  = "Ce";
        nameMap["LA"]  = "La";
        nameMap["UO2"] = "UO2";
        nameMap["CD"]  = "Cd";
        nameMap["AG"]  = "Ag";
        nameMap["BO2"] = "BO2";
        nameMap["H2O"] = "H2O";
        //nameMap["CON"] = "Concret";
        nameMap["CSI"] = "CsI";
        nameMap["CSM"] = "Cs2MoO4";
        nameMap["HE"]  = "He";
        nameMap["H2"]  = "H2";
        nameMap["O2"]  = "O2";
        nameMap["N2"]  = "N2";
        nameMap["H2O"] = "H2O";
        nameMap["H2O-VAP"] = "H2O";
        nameMap["ZR"] = "Zr";

        nameMap["LI"]  = "Li";
        nameMap["I"]   = "I";
        nameMap["TH"]  = "Th";
        nameMap["U"]   = "U";




        // from GEMS to MELCOR
        nameMap["Xe"]  = "XE";
        nameMap["Cs"]  = "CS";
        nameMap["Ba"]  = "BA";
        nameMap["I2"]  = "I2";
        nameMap["I"]   = "I2";
        nameMap["Te"]  = "TE";
        nameMap["Ru"]  = "RU";
        nameMap["Mo"]  = "MO";
        nameMap["Ce"]  = "CE";
        nameMap["La"]  = "LA";
        nameMap["UO2"] = "UO2";
        nameMap["U"]   = "U";
        nameMap["Cd"]  = "CD";
        nameMap["Ag"]  = "AG";
        nameMap["BO2"] = "BO2";
        nameMap["H2O"] = "H2O";
        //nameMap["Concret"] = "CON";
        nameMap["CsI"] = "CSI";
        nameMap["Cs2MoO4"] = "CSM";
        nameMap["He"]  = "HE";
        nameMap["H2"]  = "H2";
        nameMap["O2"]  = "O2";
        nameMap["N2"]  = "N2";
        nameMap["Zr"]  = "ZR";
        nameMap["Li"]  = "LI";
        nameMap["Th"]  = "TH";


        diff["Xe"].k0  = 2.00e5;
        diff["Cs"].k0  = 2.00e5;
        diff["Ba"].k0  = 2.95e5;
        diff["I2"].k0  = 2.00e5;
        diff["I"].k0   = 2.00e5;
        diff["Te"].k0  = 2.00e5;
        diff["Ru"].k0  = 1.62e6;
        diff["Mo"].k0  = 2.35e1;
        diff["Ce"].k0  = 2.67e8;
        diff["La"].k0  = 1.46e7;
        diff["UO2"].k0 = 1.46e7;
        diff["U"].k0   = 1.46e7;
        diff["O"].k0   = 1.46e7;
        diff["Cd"].k0  = 5.95e3;
        diff["Ag"].k0  = 5.95e3;
        diff["BO2"].k0 = 2.95e5;
        diff["B"].k0   = 2.95e5;
        diff["H2O"].k0 = 2.95e5;
        diff["H"].k0   = 2.95e5;
        //diff["CON"].k0 = 2.95e5;
        diff["CsI"].k0 = 2.00e5;
        diff["CsM"].k0 = 2.00e5;
        diff["Zr"].k0  = 0.00e0;

        diff["Xe"].Q  = 63.8;
        diff["Cs"].Q  = 63.8;
        diff["Ba"].Q  = 100.2;
        diff["I2"].Q  = 63.8;
        diff["I"].Q  = 63.8;
        diff["Te"].Q  = 63.8;
        diff["Ru"].Q  = 152.8;
        diff["Mo"].Q  = 44.1;
        diff["Ce"].Q  = 188.2;
        diff["La"].Q  = 143.1;
        diff["U"].Q   = 143.1;
        diff["O"].Q   = 143.1;
        diff["UO2"].Q = 143.1;
        diff["Cd"].Q  = 70.8;
        diff["Ag"].Q  = 70.8;
        diff["BO2"].Q = 100.2;
        diff["B"].Q   = 100.2;
        diff["H2O"].Q = 55.0;
        diff["H"].Q   = 100.0;
        //diff["CON"].Q = 100.0;
        diff["CsI"].Q = 63.8;
        diff["CsM"].Q = 63.8;
        diff["Zr"].Q  = 35.0;

        molMass["I2"]  = 126.90447 * 2;
        molMass["UO2"] = 238.0289  + 15.9994 * 2;
        molMass["BO2"] = 10.811    + 15.9994 * 2;
        molMass["H2O"] = 1.00794*2 + 15.9994;
        molMass["CsI"] = 132.90543 + 126.90447;
        molMass["Cs2MoO4"] = 132.90543*2 + 95.94 + 15.9994 * 4;
        molMass["H2"]  = 1.00794*2;
        molMass["O2"]  = 15.9994 * 2;
        molMass["N2"]  = 14.00674 * 2;
        //molMass["Concret"]  = 1.00000000;  // Not known yet
        molMass["CSM"] = 132.90543*2 + 95.94 + 15.9994 * 4;
        molMass["CSI"] = 132.90543 + 126.90447;

        molMass["H"]   = 1.00794;
        molMass["He"]  = 4.002602;
        molMass["Li"]  = 6.941;
        molMass["Be"]  = 9.012182;
        molMass["B"]   = 10.811;
        molMass["C"]   = 12.011;
        molMass["N"]   = 14.00674;
        molMass["O"]   = 15.9994;
        molMass["F"]   = 18.9984032;
        molMass["Ne"]  = 20.1797;
        molMass["Na"]  = 22.989768;
        molMass["Mg"]  = 24.305;
        molMass["Al"]  = 26.981539;
        molMass["Si"]  = 28.0855;
        molMass["P"]   = 30.973762;
        molMass["S"]   = 32.066;
        molMass["Cl"]  = 35.4527;
        molMass["Ar"]  = 39.948;
        molMass["K"]   = 39.0983;
        molMass["Ca"]  = 40.078;
        molMass["Sc"]  = 44.95591;
        molMass["Ti"]  = 47.88;
        molMass["V"]   = 50.9415;
        molMass["Cr"]  = 51.9961;
        molMass["Mn"]  = 54.93805;
        molMass["Fe"]  = 55.847;
        molMass["Co"]  = 58.9332;
        molMass["Ni"]  = 58.6934;
        molMass["Cu"]  = 63.546;
        molMass["Zn"]  = 65.39;
        molMass["Ga"]  = 69.723;
        molMass["Ge"]  = 72.61;
        molMass["As"]  = 74.92159;
        molMass["Se"]  = 78.96;
        molMass["Br"]  = 79.904;
        molMass["Kr"]  = 83.8;
        molMass["Rb"]  = 85.4678;
        molMass["Sr"]  = 87.62;
        molMass["Y"]   = 88.90585;
        molMass["Zr"]  = 91.224;
        molMass["ZR"]  = 91.224;
        molMass["Nb"]  = 92.90638;
        molMass["Mo"]  = 95.94;
        molMass["MO"]  = 95.94;
        molMass["Tc"]  = 97.9072;
        molMass["Ru"]  = 101.07;
        molMass["RU"]  = 101.07;
        molMass["Rh"]  = 102.9055;
        molMass["Pd"]  = 106.42;
        molMass["Ag"]  = 107.8682;
        molMass["AG"]  = 107.8682;
        molMass["Cd"]  = 112.411;
        molMass["CD"]  = 112.411;
        molMass["In"]  = 114.818;
        molMass["Sn"]  = 118.71;
        molMass["Sb"]  = 121.757;
        molMass["Te"]  = 127.6;
        molMass["TE"]  = 127.6;
        molMass["I"]   = 126.90447;
        molMass["Xe"]  = 131.29;
        molMass["XE"]  = 131.29;
        molMass["Cs"]  = 132.90543;
        molMass["CS"]  = 132.90543;
        molMass["Ba"]  = 137.327;
        molMass["BA"]  = 137.327;
        molMass["La"]  = 138.9055;
        molMass["LA"]  = 138.9055;
        molMass["Ce"]  = 140.115;
        molMass["CE"]  = 140.115;
        molMass["Pr"]  = 140.90765;
        molMass["Nd"]  = 144.24;
        molMass["Pm"]  = 144.9127;
        molMass["Sm"]  = 150.36;
        molMass["Eu"]  = 151.965;
        molMass["Gd"]  = 157.25;
        molMass["Tb"]  = 158.92534;
        molMass["Dy"]  = 162.5;
        molMass["Ho"]  = 164.93032;
        molMass["Er"]  = 167.26;
        molMass["Tm"]  = 168.93421;
        molMass["Yb"]  = 173.04;
        molMass["Lu"]  = 174.967;
        molMass["Hf"]  = 178.49;
        molMass["Ta"]  = 180.9479;
        molMass["W"]   = 183.84;
        molMass["Re"]  = 186.207;
        molMass["Os"]  = 190.23;
        molMass["Ir"]  = 192.22;
        molMass["Pt"]  = 195.08;
        molMass["Au"]  = 196.96654;
        molMass["Hg"]  = 200.59;
        molMass["Tl"]  = 204.3833;
        molMass["Pb"]  = 207.2;
        molMass["Bi"]  = 208.98037;
        molMass["Po"]  = 208.9824;
        molMass["At"]  = 209.9871;
        molMass["Rn"]  = 222.0176;
        molMass["Fr"]  = 223.0197;
        molMass["Ra"]  = 226.0254;
        molMass["Ac"]  = 227.0278;
        molMass["Th"]  = 232.0381;
        molMass["Pa"]  = 231.03588;
        molMass["U"]   = 238.0289;
        molMass["Np"]  = 237.0482;
        molMass["Pu"]  = 244.0642;
        molMass["Am"]  = 243.0614;
        molMass["Cm"]  = 247.0703;
        molMass["Bk"]  = 247.0703;
        molMass["Cf"]  = 251.0796;
        molMass["Es"]  = 252.083;
        molMass["Fm"]  = 257.0951;
        molMass["Md"]  = 258.0984;
        molMass["No"]  = 259.1011;
        molMass["Lr"]  = 262.1098;
        molMass["Rf"]  = 261.1089;
        molMass["Ha"]  = 262.1144;
        molMass["Sg"]  = 263.1186;
        molMass["Ns"]  = 262.1231;
        molMass["Hs"]  = 265.1306;
        molMass["Mt"]  = 266.1378;
        molMass["Unn"] = 268;
        molMass["Unu"] = 269;

    }
};

struct FilesStruct {
    std::string gems;
    std::string melcor;
    std::string restart;
    std::string melout;
    std::string meltable;
    std::string melinp;
    std::string melgen;
    std::string setup;
    std::string ptf;

    void reset() {
        gems    = "";
        melcor  = "";
        melinp  = "MEL_INPUT.inp";
        melgen  = "";
        melout  = "MELOUT_v2-0";
        restart = "MELRST_v2-0";
        setup   = "setup.inp";
        meltable= "";
        ptf     = "";
    }
};

struct Controls {
    double dTime;
    double Time;

    Controls() {
        dTime = 100.0;
        Time  = 1000.0;
    }
};

struct ValueUnit {
    std::string unut;
    double value, scale, shift;
};

struct groupStr {
    std::string group;
    int sc;
    groupStr() {
        sc = 1;
        group = "";
    }

    public:
    void Clear(){
        group = "";
        sc = 1;
    }
};

struct SetupStruct {
    std::vector<std::string> inputCols;
    std::vector<std::string> outputCols;
    int population;

    void clear() {
        inputCols.clear();
        outputCols.clear();
        population = 0;
    }
};

struct D3 {
    std::vector<double> val;

    D3(){
        val.resize(3);
    }
};

struct Properties {
    double Time;
    double TK;
    unsigned long i;
    std::string comps;
    D3      pressure;
    double  P;
    double  vPT;
    double  vPA;
    double  vPB;
    double  vPC;
    double act = 0.0;

    Properties() {
        clear();
    }

    void clear() {
        comps = "";
        P  = 0.0;
        vPT = 0.0;
        vPA = 0.0;
        vPB = 0.0;
        vPC = 0.0;
        double act = 0.0;
    }

};

struct Compound {
    double mol = 0.0;
    double kg  = 0.0;
    double act = 0.0;
    std::string name = "", phase = "s";
    std::map<std::string, double> elems;
    Properties prop;
    void clear() {
        elems.clear();
        mol  = 0.0;
        act  = 0.0;
        kg   = 0.0;
        name = 0.0;
        prop.clear();
    }
};

struct Case {
    double Time;
    double TK;
    unsigned long i;
    bool calc;
    std::vector<Compound> input;
    std::vector<Compound> compounds;
    std::map<std::string, double> elems;
    std::map<std::string, double> inelems;

    Case() {
        clear();
    }

    void clear() {
        compounds.clear();
        input.clear();
        elems.clear();
        inelems.clear();
    }
};

#endif // DATATYPES_H
