#include "reader.h"
#include <string>
#include <map>
#include "vector.h"
#include <cmath>

Reader::Reader() {

    setup.clear();
    cases.clear();
    files.reset();
    parser.Clear();
}

std::vector<std::string> Reader::split(const std::string& str){
    std::vector<std::string> result;
    std::istringstream iss(str);
    for(std::string s; iss >> s; )
        result.push_back(s);
    return result;
}

bool Reader::fileExists(const std::string& filename) {
    std::ifstream infile(filename.c_str());
    return infile.good();
}

unsigned long Reader::getDCompsCount() {
    return  node->pCSD()->nDC;
}

unsigned long Reader::getICompsCount() {
    return  node->pCSD()->nIC;
}

unsigned long Reader::getPhasesCount() {
    return  node->pCSD()->nPHb;
}

std::string Reader::getDCName(unsigned int i) {
    return node->pCSD()->DCNL[i];
}

std::string Reader::getICName(unsigned int i) {
    return node->pCSD()->ICNL[i];
}

std::string Reader::getPhaseName(int i) {
    return node->pCSD()->PHNL[i];
}

std::string Reader::getPhase(const std::string& name) {
    std::string phase = "s";
    if (name.find("(g)") != std::string::npos)
        phase = "g";
    if (name.find("(l)") != std::string::npos)
        phase = "l";
    if (name.find("(s)") != std::string::npos)
        phase = "s";
    return phase;
}

std::string Reader::cleanName(const std::string& name) {
    std::string comp = name;

    if (name.find("(g)") != std::string::npos ||
        name.find("(l)") != std::string::npos ||
        name.find("(s)") != std::string::npos)
        comp = name.substr(0, name.size() - 3);
    return comp;
}



void Reader::setAmmountOfDC(std::string name, double val) {

    std::map<std::string, int> res = parser.readFormula(name);
    std::map<std::string, int>::const_iterator it;

    for (it = res.begin(); it != res.end(); it++) {
        long int b = node->IC_name_to_xCH(it->first.c_str());
        node->Set_bIC(b, val*it->second);
    }
}

void Reader::setTemperatureC(double T) {
    node->setTemperature(T+273.15);
}

void Reader::setTemperatureK(double T) {
    node->setTemperature(T);
}

void Reader::setPressureBar(double P) {
    node->setPressure(P);
}

void Reader::setPressurePa(double P) {
    node->setPressure(P/101325);
}

void Reader::setAmmountOfIC(std::string name, double val) {
    long int b = node->IC_name_to_xCH(name.c_str());
    if (b > -1) {
        if (val >= 1.0e-10)
            node->Set_bIC(b, val);
        else
            node->Set_bIC(b, 1.0e-10);
    }
}

ValueUnit Reader::getValue(std::string& str) {

    ValueUnit res;

    std::string::iterator c;
    std::string val = "";
    std::string unit = "";


    for (c = str.begin(); c < str.end(); ++c) {

        if (isdigit(*c) || (*c) == '.' || (*c) == ',') {
            val += *c;
        } else if (isalpha(*c)) {
            unit += (*c);
        } else {
        }
    }


    res.scale = 1.0;
    res.shift = 0.0;

    if (unit.size() > 0) {
        res.unut  = unit;
        res.value = ToDouble(val);

        if (res.unut == "fs")
            res.scale = 1.0;
        if (res.unut == "ps")
            res.scale = 1.0e3;
        if (res.unut == "ns")
            res.scale = 1.0e6;

        if (res.unut == "A")
            res.scale = 1.0;
        if (res.unut == "nm")
            res.scale = 1.0e1;
        if (res.unut == "pm")
            res.scale = 1.0e-2;

        if (res.unut == "K")
            res.shift = 0.0;
        if (res.unut == "C")
            res.shift = +273.15;

        if (res.unut == "MPa")
            res.scale = 1.0;
        if (res.unut == "kPa")
            res.scale = 1.0e-3;
        if (res.unut == "GPa")
            res.scale = 1.0e3;
        if (res.unut == "Bar" || res.unut == "bar")
            res.scale = 1.0e-1;
    } else {
        res.scale = 1.0;
        res.shift = 0.0;
        res.value = ToDouble(val);
    }

    res.value = (res.value + res.shift) * res.scale;

    return res;

}

double Reader::ToDouble(const std::string& s) {
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        throw std::runtime_error("convertToDouble(\"" + s + "\")");
    return x;
}

int Reader::ToInt(const string& s){
    int r = 0;
    string::const_iterator p = s.begin();
    bool neg = false;
    if (*p == '-') {
        neg = true;
        ++p;
    }
    while (p != s.end() && *p >= '0' && *p <= '9') {
        r = (r*10) + (*p - '0');
        ++p;
    }
    if (neg) {
        r = -r;
    }
    return r;
}

void Reader::pauseMELCOR() {
    std::ofstream o("PAUSEFILE");
}

void Reader::pauseMELCOR(const std::string& path) {
    std::ofstream o(path + "PAUSEFILE");
}

/*
void Reader::getElements(MelcorCase& mel) {
    std::map<std::string, double>::const_iterator it;
    std::map<std::string, int> res;
    std::map<std::string, double> elcomp;
    std::string gem;
    double mol;

    mel.surface.elems.clear();

    for (unsigned long i  = 0; i < mel.surface.name.size(); i++) {
        gem = mel.surface.name[i];
        mol = mel.surface.mol[mel.surface.name[i]];
        elcomp = getElemCompos(gem, mol);

        for (it = elcomp.begin(); it != elcomp.end(); it++) {
            addStrToVec(melcor.surf, it->first);
            mel.surface.elems[it->first] += it->second;
        }
    }

    mel.core.elems.clear();

    for (unsigned long i = 0; i < mel.core.name.size(); i++) {
        gem = mel.core.name[i];
        mol = mel.core.mol[mel.core.name[i]];
        elcomp = getElemCompos(gem, mol);
        for (it = elcomp.begin(); it != elcomp.end(); it++) {
            addStrToVec(melcor.core, it->first);
            mel.core.elems[it->first] += it->second;
        }
    }

    mel.gas.elems.clear();

    for (unsigned long i = 0; i < mel.gas.name.size(); i++) {
        gem = mel.gas.name[i];
        mol = mel.gas.mol[mel.gas.name[i]];
        elcomp = getElemCompos(gem, mol);
        for (it = elcomp.begin(); it != elcomp.end(); it++) {
            addStrToVec(melcor.gas, it->first);
            mel.gas.elems[it->first] += it->second;
        }
    }
}
*/

void Reader::getElements(Case& cas) {
    std::map<std::string, double>::const_iterator it;
    std::map<std::string, double>::const_iterator add_it;
    std::map<std::string, int> res;
    std::map<std::string, double> elcomp;
    std::string gem;
    double mol;

    cas.inelems.clear();

    for (Compound& comp : cas.input) {
        gem = cleanName(comp.name);
        mol = comp.mol;
        comp.elems = getElemCompos(gem, mol);

        for (it = comp.elems.begin(); it != comp.elems.end(); it++) {
            cas.inelems[it->first] += it->second;
        }
    }
}

void Reader::addStrToVec(std::vector<std::string>& vec, std::string str) {
    if (std::find(vec.begin(), vec.end(), str) == vec.end())
        vec.push_back(str);
}

double Reader::getMolarMass(const std::string& name) {
    std::map<std::string, int>::const_iterator it;
    double tot;
    std::string comp = name;

    std::map<std::string, int> res = parser.readFormula(comp);
    tot = 0.0;
    for (it = res.begin(); it != res.end(); it++) {
        tot += data.molMass[it->first] * it->second;
    }
    return tot;
}

std::map<std::string, double> Reader::getElemCompos(const std::string& name, const double& mol) {
    std::map<std::string, int>::const_iterator it;
    std::map<std::string, int> res;
    std::map<std::string, double> compos;
    std::string comp = name;

    double tot, sc;
    compos.clear();

    tot = getMolarMass(comp);
    sc = mol / tot;

    res = parser.readFormula(comp);
    for (it = res.begin(); it != res.end(); it++) {
        compos[it->first] += it->second * mol;
    }

    return compos;
}

double Reader::elementChemPotential(std::string& el) {
    long int b;
    b = node->IC_name_to_xCH(el.c_str());
    return node->pCNode()->uIC[ b ];
}

double Reader::dcActivity(std::string& el) {
    long int b;
    b = node->DC_name_to_xCH(el.c_str());
    return node->pCNode()->gam[ b ];
}

std::map<std::string, double> Reader::getElemActivities(const std::string& name, const double& mol) {
    std::map<std::string, int>::const_iterator it;
    std::map<std::string, int> res;
    std::map<std::string, double> compos;
    std::string comp = name;

    return compos;
}


void Reader::setCase(Case& mel) {
    std::map<std::string, double>::const_iterator it;

    getElements(mel);

    for (unsigned int i = 0; i < getICompsCount(); i++) {
        setAmmountOfIC(getICName(i), 1.0e-8);
    }

    setTemperatureK(mel.TK);

    for (it = mel.inelems.begin(); it != mel.inelems.end(); it++) {
        if (it->second > 1.0e-8)
            setAmmountOfIC(it->first, it->second * data.sc);
    }

    /*
    if (cases.size()>1) {
        for (it = cases[1].inelems.begin(); it != cases[1].inelems.end(); it++) {
            if (it->second > 1.0e-8)
                setAmmountOfIC(it->first, it->second);
        }
    } else {
        for (it = mel.inelems.begin(); it != mel.inelems.end(); it++) {
            if (it->second > 1.0e-8)
                setAmmountOfIC(it->first, it->second);
        }
    }
    */
    /*
    for (it = mel.inelems.begin(); it != mel.inelems.end(); it++) {
        if (it->second > 1.0e-8)
            setAmmountOfIC(it->first, it->second);
    }
    */
}

void Reader::printMELCOR(Case& mel, std::string message) {
    std::map<std::string, double>::const_iterator it;
    std::map<std::string, double> elems;
    std::map<std::string, double> tot;
    std::stringstream stream;

    for (it = mel.elems.begin(); it != mel.elems.end(); it++) {
        elems[it->first] += mel.elems[it->first];
    }

    std::cout << message << "   " << mel.TK << "   "  << mel.i << std::endl;
    for (it = mel.elems.begin(); it != mel.elems.end(); it++) {
        tot[it->first] = mel.elems[it->first];
    }

        stream.str("");
        stream.clear();
        stream.fill(' ');
        stream.width(15);
        stream.precision(8);
        stream << scientific << tot["Cs"] << "\t" << scientific << elems["Cs"];
        std::cout << "Cs" << "\t" << scientific << stream.str() << std::endl;
}

/*
void Reader::saveDF() {
    std::stringstream stream;
    std::fstream dfFile;

    //std::cout << melcor[0].surface.name.size() << std::endl;

    dfFile.open("result.df", std::ios::out);
    std::map<std::string, std::vector<double>>::const_iterator iter;

    stream.str("");
    stream.clear();
    stream.fill(' ');

    stream << "t" << ";" << "tk" << ";";


    // add data on core composition
    for (std::string st : melcor.core){
        stream << st << "(core)" << ";" ;
    }

    // add data on surface composition
    for (std::string st : melcor.surf){
        stream << st << "(surf)" << ";" ;
    }

    // add data on surface composition
    for (std::string st : melcor.gas){
        stream << st << "(gas)" << ";" ;
    }

    for (Compound st : gems.cases[0].compounds){
        stream << st.name << ";" ;
    }

    dfFile << stream.str() << std::endl;

    for ( int i = 0; i < melcor.cases.size(); i++ ) {

        stream.str("");
        stream.clear();
        stream.fill(' ');
        stream << melcor.cases[i].Time << ";";
        stream << melcor.cases[i].TK << ";";



        for (Compound st : gems.cases[i].compounds){
            stream << st.mol << ";" ;
        }

        dfFile << stream.str() << std::endl;
    }

    dfFile.close();
}
*/

bool Reader::findNextLine(std::fstream& file, std::string& find, std::string& line) {
    std::vector<std::string> res;
    std::vector<std::string> comp;

    comp = split(find);

    std::getline(file, line);
    res   = split(line);

    bool is_equal;

    is_equal = false;
    if(comp.size() == res.size())
        is_equal = std::equal(comp.begin(), comp.end(), res.begin());
    if (is_equal) return true;

    while (!file.eof()) {
        std::getline(file, line);
        res   = split(line);
        is_equal = false;
        if(comp.size() == res.size())
            is_equal = std::equal(comp.begin(), comp.end(), res.begin());
        if (is_equal) return true;
    }
    return false;
}

bool Reader::findNextLinePart(std::fstream& file, std::string& find, std::string& line) {
    std::vector<std::string> res;
    std::vector<std::string> comp;

    std::getline(file, line);

    if (isEqual(find, line)) return true;

    while (!file.eof()) {
        std::getline(file, line);
        if (isEqual(find, line)) return true;
    }
    return false;
}

bool Reader::isEqual(std::string& find, std::string& line) {
    std::vector<std::string> res;
    std::vector<std::string> comp;
    bool is_equal;

    comp = split(find);
    res   = split(line);

    is_equal = false;
    if(comp.size() <= res.size())
        is_equal = std::equal(comp.begin(), comp.end(), res.begin());
    return is_equal;
}

/*
void Reader::diffuse() {
    double dval, dt, k0, Q, f, time, TimeStep;
    double R = 1.987e-3;
    dt = 0.01;

    //Modified CORSOR-Booth
    double pi, D, D_dim, Q_CB, a;
    pi = 3.142;
    a = 6.0E-6;    //radius of the fuel grain, MELCOR default = 6.0E-6
    double ra2 = 1.0 / (a*a);
    D0 = 1.0E-6		// MELCOR default
    Q_CB = 3.814E5 	//MELCOR default

    for ( DataMelcor& mel : melcor ) {

        if (mel.i > 0) {

            for (unsigned long i = 0; i < mel.core.kg.size(); i++) {
                mel.core.kg[i]    = melcor[mel.i-1].core.kg[i];
                mel.surface.kg[i] = melcor[mel.i-1].surface.kg[i];
            }

            for (unsigned long i = 0; i < mel.core.kg.size(); i++) {
                TimeStep = mel.Time - melcor[mel.i-1].Time;

                // CORSOR-M
                k0 = data.k0[mel.core.name[i]] / 60.0;
                Q  = data.Q[mel.core.name[i]];
                f  = k0 * exp( - Q / ( R * mel.TK ));

                //Modifies CORSOR-Booth
                D0 = data.D0[mel.core.name[i]];
                Q_CB  = data.Q_CB[mel.core.name[i]];
                S_k  = data.S_k[mel.core.name[i]];
                D  = D0 * exp( - Q / ( R * mel.TK ));
                D_dim = D * ra2;


                time = 0.0;

                while (time < TimeStep) {

                    //CORSOR-M
                    dval = f * dt * (mel.core.kg[i] - mel.surface.kg[i]);
                    mel.core.kg[i]    -= dval;
                    mel.surface.kg[i] += dval;

                    //Modified CORSOR-Booth
                    if (D_dim * time < 1/pow(pi,2)) {
                        f = 6 * sqrt(D_dim*time/pi) - 3*D_dim*time;
                    }
                    else {
                        f = 1 - 6/pow(pi,2)*exp(-pow(pi,2)*D_dim*time);
                    }
                    Rel_Cs = (f[k] - f[k-1]) / dt; // not work Ja
                    Rel_k = S_k * Rel_Cs;
                    dval = Rel_k * dt * (mel.core.kg[i] - mel.surface.kg[i]);
                    mel.core.kg[i]    -= dval;
                    mel.surface.kg[i] += dval;

                    time += dt;

                }
            }
        } else {
            for (unsigned long i = 0; i < mel.core.kg.size(); i++) {
                mel.surface.kg[i] = 0.0;
            }
        }
    }
}
*/


/*
void Reader::diffuse(MelcorCase& mel) {
    std::map<std::string, double>::iterator it;
    double dval, dt, dtk, k0, Q, f, time, TimeStep;
    double R = 1.987e-3;
    double steps, TK;
    dt = 0.1;

    if (mel.i > 0) {

        steps = (mel.Time - melcor.cases[mel.i-1].Time) / dt;
        dtk   = (mel.TK   - melcor.cases[mel.i-1].TK)     / steps;
        TimeStep = mel.Time - melcor.cases[mel.i-1].Time;

        mel.core.elems    = melcor.cases[mel.i-1].core.elems;
        mel.surface.elems = melcor.cases[mel.i-1].surface.elems;

        //getElements(mel);
        //printMELCOR(mel, " ---------- BEFORE DIFFUSION -----------");

        for (it = mel.core.elems.begin(); it != mel.core.elems.end(); it++) {
            k0 = data.diff[it->first].k0 / 60.0;
            Q  = data.diff[it->first].Q;

            time  = 0.0;
            TK    = mel.TK;

            while (time < TimeStep) {
                f  = k0 * exp( - Q / ( R * TK ));
                dval = f * dt * ( it->second - mel.surface.elems[it->first] );
                //dval = f * dt * mel.core.mol[i];
                it->second  -= dval;
                mel.surface.elems[it->first] += dval;
                time += dt;
                TK   += dtk;
            }
        }
    }
}


void Reader::diffuse(MelcorCase& mel, CompData& from, CompData& to) {
    std::map<std::string, double>::iterator it;
    std::map<std::string, double> cfrom;
    std::map<std::string, double> cto;
    double dval, dt, dtk, k0, Q, f, time, TimeStep;
    double R = 1.987e-3;
    double steps, TK, totfrom, totto;
    std::string comp;
    dt = 0.1;

    if (mel.i > 0) {

        steps    = (mel.Time - melcor.cases[mel.i-1].Time) / dt;
        dtk      = (mel.TK   - melcor.cases[mel.i-1].TK)   / steps;
        //TimeStep = (mel.Time - melcor.cases[mel.i-1].Time);

        totfrom = 0;
        for (it = from.mol.begin(); it != from.mol.end(); it++)
            if (!std::isnan(it->second))
                totfrom += it->second;
            else
                totfrom += 0.0;
        for (it = from.mol.begin(); it != from.mol.end(); it++)
            if (!std::isnan(it->second))
                cfrom[it->first] = it->second / totfrom;
            else
                cfrom[it->first] = 0.0;

        totto = 0.0;
        for (it = to.mol.begin(); it != to.mol.end(); it++)
            totto += it->second;
        for (it = to.mol.begin(); it != to.mol.end(); it++)
            cto[it->first] = it->second / totto;
        //from.elems   = melcor.cases[mel.i-1].core.elems;
        //to.elems     = melcor.cases[mel.i-1].surface.elems;

        //printMELCOR(mel, " ---------- BEFORE DIFFUSION -----------");

        for (it = from.mol.begin(); it != from.mol.end(); it++) {
            comp = cleanName(it->first);
            k0 = data.diff[comp].k0 / 60.0;
            Q  = data.diff[comp].Q;

            time  = 0.0;
            TK    = mel.TK;
            double& fr = it->second;
            double& tt = to.mol[it->first];
            double& cf = cfrom[it->first];
            double& ct = cto[it->first];

            //while (time < TimeStep) {
            while (time < control.dTime) {
                f  = k0 * exp( - Q / ( R * TK ));
                dval = f * dt * ( fr - tt );
                fr  -= dval;
                tt += dval;
                cf  = fr / totfrom;
                ct  = tt / totto;
                time += dt;
                TK   += dtk;
            }

        }
    }
}
*/

double Reader::getTboil(std::string comp) {
    long int l = node->DC_name_to_xCH((comp+"(l)").c_str()); // code for the liquid phase of the compound
    long int g = node->DC_name_to_xCH((comp+"(g)").c_str()); // code for the gas phase of the compound
    double T1=273.15, T2 = 273.15+80.0, P = 100000.0, dG1, dG2, dT = 80.0, epsilon = 1.0e-2, err = 1.0;
    bool norm = false; // in this case the DC_GO will be in J/mol units

    if (l < 0 || g < 0) return 1.0;

    dG1 = node->DC_G0(g, P, T1, norm) - node->DC_G0(l, P, T1, norm);
    dG2 = node->DC_G0(g, P, T2, norm) - node->DC_G0(l, P, T2, norm);

    while (err > epsilon) {
        if (dG1 > 0.0 & dG2 > 0.0 ) {
            T1 = T2; 
            T2 +=dT;
        }
        if (dG1 < 0.0 & dG2 < 0.0 ) {
            T2 = T1; 
            T1 -= dT;
        }
        if (dG1 > 0.0 & dG2 < 0.0 ) {
            dT *= 0.5;
            T2 = T1; 
            T1 = T1 - dT;
        }

        dG1 = node->DC_G0(g, P, T1, norm) - node->DC_G0(l, P, T1, norm);
        dG2 = node->DC_G0(g, P, T2, norm) - node->DC_G0(l, P, T2, norm);
        err = abs(T1 - T2);

        if (T1 >= 2273.0) return 2273.0;
    }
    return (T1 + T2) * 0.5;
}

double Reader::getTmelt(std::string comp) {
    long int s = node->DC_name_to_xCH((comp+"").c_str()); // code for the liquid phase of the compound
    long int l = node->DC_name_to_xCH((comp+"(l)").c_str()); // code for the gas phase of the compound
    double T1=273.15, T2 = 273.15+80.0, P = 100000.0, dG1, dG2, dT = 80.0, epsilon = 1.0e-2, err = 1.0;
    bool norm = false; // in this case the DC_GO will be in J/mol units

    if (l < 0 || s < 0) return 1.0;

    dG1 = node->DC_G0(l, P, T1, norm) - node->DC_G0(s, P, T1, norm);
    dG2 = node->DC_G0(l, P, T2, norm) - node->DC_G0(s, P, T2, norm);

    while (err > epsilon) {
        if (dG1 > 0.0 & dG2 > 0.0 ) {
            T1 = T2;
            T2 +=dT;
        }
        if (dG1 < 0.0 & dG2 < 0.0 ) {
            T2 = T1;
            T1 -= dT;
        }
        if (dG1 > 0.0 & dG2 < 0.0 ) {
            dT *= 0.5;
            T2 = T1;
            T1 = T1 - dT;
        }

        dG1 = node->DC_G0(l, P, T1, norm) - node->DC_G0(s, P, T1, norm);
        dG2 = node->DC_G0(l, P, T2, norm) - node->DC_G0(s, P, T2, norm);
        err = abs(T1 - T2);

        if (T1 >= 2273.0) return 2273.0;
    }
    return (T1 + T2) * 0.5;
}

bool Reader::contains(std::string str, std::string f) {
    return str.find(f) != std::string::npos;
}

int Reader::getIndex(std::vector<Compound>& vec, std::string& str, std::string& phase) {
    int pos = -1;
    int i = 0;
    for (auto comp : vec) {
        if (comp.name == str && comp.phase == phase) {
            pos = i;
            break;
        }
        i++;
    }
    return pos;
}


std::vector<double> Reader::CCPressureCoeffs(std::string& name) {
    std::string comp;
    double Tb, DHvap;
    int g, l;
    std::vector<double> res(3);

    if (contains(name, "(g)") || contains(name, "(l)") || contains(name, "(s)")) {
        unsigned int pos = name.size() - 3;
        comp = name.erase(pos, 3);
    } else {
        comp = name;
    }

    Tb     = getTboil(comp);
    g = node->DC_name_to_xCH((comp+"(g)").c_str()); // code for the liquid phase of the compound
    l = node->DC_name_to_xCH((comp+"(l)").c_str()); // code for the gas phase of the compound

    DHvap  = node->DC_H0(g, 1.0e5, Tb) - node->DC_H0(l, 1.0e5, Tb);

    res[0] = 0.43429 *  DHvap / 8.3145;
    res[1] = 2.88081 + 0.43429 * DHvap / (8.3145 * Tb); // 2.88081 is a scaling factor from bar to mmHg (in log10 units) log10(760)
    res[2] = 0.0;

    return res;
}

/*
/// calculated the vapour pressure based on the Clapeyron-Clausius expression from the Hgas-Hliq
void Reader::getPressure(Case& gems) {
    std::map<std::string, int>::const_iterator it;
    std::string name, comp="LiF", phase="g";
    double mol, lga, Tb, DHvap, t;
    long int g, l, ind;
    std::map<std::string, double> elems, melc, condens;
    std::map<std::string, double>::iterator elit;
    std::map<std::string, int> res;
    std::vector<double> y, b(3), x;

        
    for (unsigned int i = 0; i < gems.compounds.size(); i++) {
        b = CCPressureCoeffs(gems.compounds[i].name);
        gems.compounds[i].prop.vPA = b[0];
        gems.compounds[i].prop.vPB = b[1];
        gems.compounds[i].prop.vPC = b[2];
        gems.compounds[i].prop.vPT = gems.compounds[i].prop.vPA / ( 10.0 + gems.compounds[i].prop.vPB ); // Temperature at which log10(P) is equat -10.0 (or pressure is 1.0e-10 mmHg);
    }
}
*/

/// calculated the vapour pressure based on the Clapeyron-Clausius expression from the Hgas-Hliq
void Reader::getPressure(Case& gems) {
    std::map<std::string, int>::const_iterator it;
    std::string name, comp="LiF", phase="g";
    double mol, lga, Tb, DHvap, t;
    long int g, l, ind;
    std::map<std::string, double> elems, melc, condens;
    std::map<std::string, double>::iterator elit;
    std::map<std::string, int> res;
    std::vector<double> y, b(3), x;

    if (cases.size() < 2) {
        
        for (unsigned int i = 0; i < gems.compounds.size(); i++) {
            b = CCPressureCoeffs(gems.compounds[i].name);
            gems.compounds[i].prop.vPA = b[0];
            gems.compounds[i].prop.vPB = b[1];
            gems.compounds[i].prop.vPC = b[2];
            gems.compounds[i].prop.vPT = gems.compounds[i].prop.vPA / ( 10.0 + gems.compounds[i].prop.vPB ); // Temperature at which log10(P) is equat -10.0 (or pressure is 1.0e-10 mmHg);
        }

    } else {
        
        for (unsigned int i = 0; i < gems.compounds.size(); i++) {
            x.resize(0); 
            y.resize(0); 

            for (unsigned int n = (cases.size()-1); n > 0; n--) {

                ind = getIndex(cases[n].compounds, gems.compounds[i].name, gems.compounds[i].phase);
                if (ind > -1) {
                    //if (cases[n].compounds[ind].prop.P > 1.0e-20 && cases[n].compounds[ind].prop.P <= 760){
                    y.push_back(log10(cases[n].compounds[ind].prop.P));
                    x.push_back(-1.0/cases[n].TK);
                    //}
                }
            }

            if (y.size() > 2) {
                double n = y.size();

                double avgX = std::accumulate(x.begin(), x.end(), 0.0) / n;
                double avgY = std::accumulate(y.begin(), y.end(), 0.0) / n;

                double numerator = 0.0;
                double denominator = 0.0;

                for(int i=0; i<n; ++i){
                    numerator += (x[i] - avgX) * (y[i] - avgY);
                    denominator += (x[i] - avgX) * (x[i] - avgX);
                }

                b[0] = numerator / denominator;
                b[1] = avgY - b[0] * avgX;
                b[2] = 0.0;

            } else {
                b = CCPressureCoeffs(gems.compounds[i].name);
            }

            gems.compounds[i].prop.vPA = b[0];
            gems.compounds[i].prop.vPB = b[1];
            gems.compounds[i].prop.vPC = b[2];
            gems.compounds[i].prop.vPT = gems.compounds[i].prop.vPA / ( 10.0 + gems.compounds[i].prop.vPB ); // Temperature at which log10(P) is equat -10.0 (or pressure is 1.0e-10 mmHg);
        }
    }
}


void Reader::rePopulateData() {
    Case mel0 = cases[0];
    std::vector<Case> buff;
    double dt, time, dtk, tk;
    unsigned long up, down, i;

    if (setup.population > 0) {
        buff.resize(setup.population);

        mel0 = cases[0];

        time = 0.0;
        tk = cases[0].TK;

        down = 0;
        up   = 1;
        dt   = ( cases[cases.size() - 1].Time - cases[0].Time ) / setup.population;
        dtk  = (cases[up].TK - cases[down].TK) / (cases[up].Time - cases[down].Time);

        buff[0] = cases[0];

        buff[0].compounds = mel0.compounds;

        i = 0;
        for ( Case& b : buff ) {
            b       = cases[down];
            b.Time  = time;
            b.TK    = tk;
            b.i     = i;

            if (time >= cases[down].Time && time < cases[up].Time) {
            } else {
                for (unsigned long i = 0; i < cases.size()-1; i++) {
                    if (time >= cases[i].Time && time < cases[i+1].Time) {
                        down = i;
                        up   = i+1;
                        dtk  = (cases[up].TK - cases[down].TK) / (cases[up].Time - cases[down].Time);
                        break;
                    }
                }
            }

            tk   += dtk * dt;
            time += dt;
            i    += 1;
        }

        cases = buff;
    }
}

std::vector<double> Reader::calcCoeffs(std::vector<std::vector<double>> x, std::vector<double> y) {
    Eigen::MatrixXd xx(x.size(),3);
    Eigen::VectorXd yy(y.size());

    Eigen::VectorXd b(3);

    yy = Eigen::VectorXd::Map(&y[0], y.size());

    for (int i = 0; i < x.size(); i++) {
        xx.row(i) = Eigen::VectorXd::Map(&x[i][0],x[i].size());
    }

    b = xx.ldlt().solve(yy); //xx.transpose() * xx).ldlt(); // .ldlt().solve(xx.transpose() * yy);

    vector<double> v;

    for (int i = 0; i < b.size(); i++) {
        v.push_back(b[i]);
    }
    return v;
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::string formatValue(double value, int prec = 0, int width = 10) {
    std::stringstream stream;
    std::string res;

    stream.clear();
    stream.fill(' ');
    if (prec > 0)
        stream.precision(prec);
    if (std::abs(value) > 0.0) {
        if (std::abs(value) > -1.0e5 && std::abs(value) < 1.0e5) {
            stream << std::right << std::fixed << std::setw(width) << value;
        } else {
            stream << std::right << std::scientific << std::setw(width) << value;
        }
    } else {
        stream << std::right << std::fixed << std::setw(width) << 0.0;
    }
    res = stream.str();
    replace(res, "e+", "E");
    replace(res, "e-", "E-");

    int size = res.length();
    return std::string(width-size, ' ') + res + "  ";
};

std::string formatName(std::string name, int width) {
    std::stringstream stream;
    std::string res, sp;


    stream.clear();
    stream.fill(' ');
    stream << "'" << name << "'";
    res = stream.str();
    int size = res.length();
    sp = std::string(width-size, ' ');
    return res + sp;
};

void Reader::writeMelcorInput(Case& gems) {
    std::ofstream file;
    std::string comp, phase = "g";
    std::stringstream stream, final, species;

    file.open (files.melinp);

    stream.clear();
    species.clear();
    final.clear();
    //stream << "MEL_OUTPUTFILE 'MC_V-1_output.out' \n";
    stream << "MEL_OUTPUTFILE '" + files.melout  + "' NCYCLE -1 \n";
    //stream << "MEL_RESTARTFILE 'V-1_restart.rst' NCYCLE -1 \n";
    stream << "MEL_RESTARTFILE '" + files.restart + "' \n";
    stream << "PROGRAM MELCOR \n";
    stream << "EXEC_INPUT \n";
    stream << "EXEC_TITLE 'VERDON-1' \n";
    stream << "EXEC_CPULEFT 20.0 \n";
    stream << "EXEC_CPULIM 1.0E5 \n";
    stream << "EXEC_CYMESF 800 800 \n";
    stream << "EXEC_TIME 5 \n";
    stream << "\t \t1 -10000.0 5.0E-3 1.0E-6 10000.0 100.0 10000.0 1.0E10 \n";
    stream << "\t \t2  0.0     0.5    1.0E-6  2000.0  20.0 1000.0  1.0E10 \n";
    stream << "\t \t3  10000.0 0.5    1.0E-6  1500.0  20.0 1500.0  1.0E10 \n";
    stream << "\t \t4  11500.0 0.5    1.0E-6  1000.0  10.0 1000.0  1.0E10 \n";
    stream << "\t \t5  15000.0 0.5    1.0E-6   500.0  5.0   500.0  1.0E10 \n";
    stream << "EXEC_TEND " << formatValue( gems.Time, 2, 12) << "\n";  //30000.0 \n";


    stream << std::endl;
    stream << "CF_INPUT \n";
    stream << "! \n";
    int count, ind;

    for (DataInfo::diPair pair: data.pressureMap) {
        comp  = data.nameMap[pair.str];
        count = pair.i;
        ind   = getIndex(gems.compounds, comp, phase);

        if (ind > -1) {
            stream << std::left << "    " << std::setw(8) << "CF_ID"  << formatName(pair.str+"_A", 10);
            stream << std::left << std::setw(6) << count + 1 << std::setw(10) << "EQUALS" << std::endl;
            stream << std::left << "    " << std::setw(8) << "CF_SAI" << std::setw(8) << "0.0";
            stream << std::left << formatValue( gems.compounds[ind].prop.vPA , 1, 14) << std::endl;
            stream << "! \n";

            stream << std::left << "    " << std::setw(8) << "CF_ID"  << formatName(pair.str+"_B", 10);
            stream << std::left << std::setw(6) << count + 2 << std::setw(10) << "EQUALS" << std::endl;
            stream << std::left << "    " << std::setw(8) << "CF_SAI" << std::setw(8) << "0.0";
            stream << std::left << formatValue( gems.compounds[ind].prop.vPB , 1, 14) << std::endl;
            stream << "! \n";

            stream << std::left << "    " << std::setw(8) << "CF_ID"  << formatName(pair.str+"_C", 10);
            stream << std::left << std::setw(6) << count + 3 << std::setw(10) << "EQUALS" << std::endl;
            stream << std::left << "    " << std::setw(8) << "CF_SAI" << std::setw(8) << "0.0";
            stream << std::left << formatValue( gems.compounds[ind].prop.vPC , 1, 14) << std::endl;
            stream << "! \n";
        }
    }
    stream << "END PROGRAM MELCOR";

    /*
    stream << std::endl;
    stream << "CF_INPUT \n";
    stream << "! \n";
    int count = 700;
    for (unsigned int i = 0; i < mel.compounds.size(); i++) {
        if (data.nameMap[mel.compounds[i].name] != "") {
            stream << std::left << "    " << std::setw(8) << "CF_ID"  << formatName(data.nameMap[mel.compounds[i].name]+"_A", 10);
            stream << std::left << std::setw(6) << count + 1 << std::setw(10) << "EQUALS" << std::endl;
            stream << std::left << "    " << std::setw(8) << "CF_SAI" << std::setw(8) << "0.0";
            stream << std::left << formatValue( mel.measurements[i].vPA , 1, 14) << std::endl;
            stream << "! \n";

            stream << std::left << "    " << std::setw(8) << "CF_ID"  << formatName(data.nameMap[mel.compounds[i].name]+"_B", 10);
            stream << std::left << std::setw(6) << count + 2 << std::setw(10) << "EQUALS" << std::endl;
            stream << std::left << "    " << std::setw(8) << "CF_SAI" << std::setw(8) << "0.0";
            stream << std::left << formatValue( mel.measurements[i].vPB , 1, 14) << std::endl;
            stream << "! \n";

            stream << std::left << "    " << std::setw(8) << "CF_ID"  << formatName(data.nameMap[mel.compounds[i].name]+"_C", 10);
            stream << std::left << std::setw(6) << count + 3 << std::setw(10) << "EQUALS" << std::endl;
            stream << std::left << "    " << std::setw(8) << "CF_SAI" << std::setw(8) << "0.0";
            stream << std::left << formatValue( mel.measurements[i].vPC , 1, 14) << std::endl;
            stream << "! \n";
            count+=10;
        }
    }
    stream << "END PROGRAM MELCOR";
    */

    /*stream << "COR_INPUT \n";
    stream << "!COR_TST 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 this disables the oxiation model \n";
    stream << "RN1_INPUT \n";

    int count = 0;

    for (unsigned int i = 0; i < mel.compounds.size(); i++) {
        if (data.nameMap[mel.compounds[i].name] != "") {
            species << std::left << "        " << std::setw(5) << count*6+1 << std::setw(8) << "7110" << formatName(data.nameMap[mel.compounds[i].name], 10);
            species << formatValue( mel.measurements[i].vPT , 1, 12);
            species << "1 1 ! temperature T (< T vapor pressure is 0) \n";

            species << std::left << "        " << std::setw(5) << count*6+2 << std::setw(8) << "7110" << formatName(data.nameMap[mel.compounds[i].name], 10);
            species << formatValue( mel.measurements[i].vPA , 4, 12);
            species << "1 2 ! constant A \n";

            species << std::left << "        " << std::setw(5) << count*6+3 << std::setw(8) << "7110" << formatName(data.nameMap[mel.compounds[i].name], 10);
            species << formatValue( mel.measurements[i].vPB , 4, 12);
            species << "1 3 ! constant B \n";

            species << std::left << "        " << std::setw(5) << count*6+4 << std::setw(8) << "7110" << formatName(data.nameMap[mel.compounds[i].name], 10);
            species << formatValue( mel.measurements[i].vPC , 4, 12);
            species << "1 4 ! constant C \n";

            species << std::left << "        " << std::setw(5) << count*6+5 << std::setw(8) << "7110" << formatName(data.nameMap[mel.compounds[i].name], 10);
            species << formatValue( -1.0 , 1, 12);
            species << "2 1 ! constant C \n";

            species << std::left << "        " << std::setw(5) << count*6+6 << std::setw(8) << "7110" << formatName(data.nameMap[mel.compounds[i].name], 10);
            species << formatValue( -1.0 , 1, 12);
            species << "3 1 ! constant C \n";
            count++;
        }
    }

    stream << "RN1_CSC " << count*6 << "\n";
    stream << "!here we add the vapor pressure info \n";
    stream << species.str();
    stream << "END PROGRAM MELCOR";
    */

    file << stream.str();
    file.close();
}

void Reader::readMelcor( ) {
    if (files.meltable != "")
        readMelcorTable(files.meltable);
    else
        readMelcorOutput(files.melout);
}

void Reader::readMelcorTable(const std::string& file ) {
    using std::getline;
    std::map<std::string, double>::const_iterator it;
    double time = 0.0, T = 1500;
    double mol, total, act, kg, molmass;
    int i = 0, cas = 0;
    //ValueUnit ValUn;

    std::fstream melcorFile;
    std::string line, buff, phase, comp;
    std::vector<std::string> res;
    std::map<std::string, double> elcompos;

    melcorFile.open(file, std::ios::in);

    if (!melcorFile.is_open()) {
        std::cout << "The melcor file [ " << file <<  " ] is not fould " << std::endl;
        return;
    }

    /// Read the last line in the file
    std::getline(melcorFile, buff);
    if (buff.size() > 5) line = buff;
    while (!melcorFile.eof()) {
        if (buff.size() > 5) line = buff;
        std::getline(melcorFile, buff);
    }

    std::cout << line << std::endl;

    Case mel;
    res   = split(line);

    if (res.size() < setup.inputCols.size()) {
        std::cout << "Number of columns names in is not equal to the number of columns the table file" << std::endl;
        std::exit(0);
    }

    for (i = 0; i < res.size(); i++) {
        if (setup.inputCols[i] == "Time" || setup.inputCols[i] == "time")
            mel.Time = ToDouble(res[i]);

        else if (setup.inputCols[i] == "none" || setup.inputCols[i] == "out" || setup.inputCols[i] == "-")
            continue;
            
        else if (setup.inputCols[i] == "T" || setup.inputCols[i] == "TK" || setup.inputCols[i] == "Temp")
            mel.TK   = ToDouble(res[i]);

        else {
            Compound com;
            phase    = getPhase(setup.inputCols[i]);
            comp     = cleanName(setup.inputCols[i]);
            molmass  = getMolarMass(comp);
            kg       = ToDouble(res[i]);
            mol      = molmass * kg;

            com.name  = comp;
            com.phase = phase;
            com.kg    = kg;
            com.mol   = mol;
            mel.input.push_back(com);
        }
    }

    for (it = data.add.begin(); it != data.add.end(); it++) {
        Compound com;
        phase    = getPhase(it->first);
        comp     = cleanName(it->first);
        molmass  = getMolarMass(comp);
        mol      = it->second;
        kg       = mol / molmass;

        com.name  = comp;
        com.phase = phase;
        com.kg    = kg;
        com.mol   = mol;
        mel.input.push_back(com);
    }

    mel.i = cases.size();

    cases.push_back(mel);

    melcorFile.close();
    std::cout << "Done loading the melcor table file ... " << std::endl;
}

void Reader::readMelcorOutput(const std::string& file ) {
    using std::getline;
    double time = 0.0, T = 1500;
    int i = 0;
    //ValueUnit ValUn;

    std::fstream melcorFile;
    std::string line, gem;
    std::string f;
    std::vector<std::string> res;

    melcorFile.open(file, std::ios::in);

    if (!melcorFile.is_open()) {
        std::cout << "The melcor file [ " << file <<  " ] is not fould " << std::endl;
        return;
    }

    std::getline(melcorFile, line);

    while (!melcorFile.eof()) {
        Case mel;
        // Find the line that is one above the time
        f = "1* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ";
        if (findNextLine(melcorFile, f, line)) {
            std::getline(melcorFile, line);
            res   = split(line);
            time = ToDouble(res[1]);
        } else {
            break;
        }

        /*
        // Find the line that indicates the beginning of the table
        f = "VOLUME MATERIAL MASS DENSITY E PPART ECV";
        if (findNextLine(melcorFile, f, line)) {
            f = "Level_3";
            if (findNextLinePart(melcorFile, f, line)) {
                f = "Level_4";
                while (!isEqual(f, line)) {
                    res   = split(line);
                    if (res[0] == "Level_3") {
                        gem = data.nameMap[res[2]];

                        addStrToVec(melcor.gas, gem);
                        addStrToVec(melcor.core, gem);
                        addStrToVec(melcor.surf, gem);

                        mel.name.push_back(gem);
                        mel.kg.push_back(data.scale[res[2]].val * ToDouble(res[3]));
                        mel.mol[gem]= (data.scale[res[2]].val * ToDouble(res[3])/data.molMass[gem]*1000);
                        std::getline(melcorFile, line);
                    } else {
                        gem = data.nameMap[res[1]];
                        addStrToVec(melcor.gas, gem);
                        addStrToVec(melcor.core, gem);
                        addStrToVec(melcor.surf, gem);
                        mel.name.push_back(gem);
                        mel.kg.push_back(data.scale[res[1]].val * ToDouble(res[2]));
                        mel.mol[gem]= (data.scale[res[1]].val * ToDouble(res[2])/data.molMass[gem]*1000);
                        std::getline(melcorFile, line);
                    }
                }
            } else {
                break;
            }
        } else {
            break;
        }
        */

        // Find the line that indicates the beginning of the temperature table
        f = "EDIT OF CORE CELL COMPONENT TEMPERATURES (K)";
        if (findNextLine(melcorFile, f, line)) {
            f = "1  3";
            if (findNextLinePart(melcorFile, f, line)) {
                res   = split(line);
                T = ToDouble(res[2]);
            }
        } else {
            break;
        }

        // Find the line that indicates the beginning of the table
        //f = "RELEASED RADIOACTIVE RADIONUCLIDE MASS AUDIT - MASSES IN KG";
        f = "RADIOACTIVE RADIONUCLIDE MASS DISTRIBUTION IN KG";
        if (findNextLine(melcorFile, f, line)) {
            for (int num_lines = 0; num_lines < 3; ++num_lines){std::getline(melcorFile, line);}

            std::getline(melcorFile, line);
            res   = split(line);
            //for (int num_lines = 0; num_lines < 17; ++num_lines){
            while (res.size() > 0) {
                Compound comp;
                // check if the key (name) exists in the map before proceeding
                if (data.nameMap.find(res[0]) != data.nameMap.end()) {
                    gem = data.nameMap[res[0]];

                    comp.name  = gem;
                    comp.phase = "s";
                    comp.kg    = ToDouble(res[1]);
                    comp.mol   = (ToDouble(res[1])/data.molMass[res[0]]*1000);
                    mel.input.push_back(comp);
                }
                std::getline(melcorFile, line);
                res   = split(line);
            }

            mel.Time = time;
            mel.TK   = T;
            mel.i    = i;
            //getElements(mel);

            cases.push_back(mel);
            i++;
        }
        std::getline(melcorFile, line);
    }

    melcorFile.close();

    std::cout << "Done loading the melcor output file ... " << std::endl;
}

void Reader::clearFiles(std::string& f) {
    std::fstream file;
    file.open(f, std::ios::out);
    file.close();
}
