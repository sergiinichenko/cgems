#include "cgems.h"
#include <iomanip>

cGems::cGems() {
}

cGems::cGems(string file) {
    cGems();
    initiateGems(file);
}

void cGems::initiateGems(std::string file) {

    node  = new TNode();
    // get the path to the file
    unsigned pos = file.find_last_of("/");
    path = file.substr(0, pos+1);

    // copy the input file to the char[] variable
    char input_system_file_list_name[256];
    int i = 0;
    for (char c : file) {
       input_system_file_list_name[i] = c;
       i++;
    }

    // (1) Initialization of GEMS3K internal data by reading  files whose names are given in the input_system_file_list_name
    int res = node->GEM_init( input_system_file_list_name );
    std::cout << "Initialization results : " << res << std::endl;
    //node->GEM_init( &data.inputPM,  &data.inputDHC, &data.inputBR);
    // Getting direct access to work node DATABR structure which exchanges the
    // data with GEM IPM3 (already filled out by reading the DBR input file)
    dBR = node->pCNode();
}


// Input file with the control settings. Has to be read after the main input files.
void cGems::controlFile(std::string cfile) {

}

void cGems::showSelected() {
    std::stringstream stream;

    stream.str("");
    stream.clear();
    stream.fill(' ');

    for (std::string st : setup.outputCols) {
        stream.width(10);
        stream << data.nameMap[st];
    }

    std::cout << stream.str() << std::endl;
    stream.str("");
    stream.clear();
    stream.fill(' ');

    for (std::string st : setup.outputCols) {
        long int b = node->DC_name_to_xCH(data.nameMap[st].c_str());
        if (node->Get_nDC(b) > 1.0e-20) {
            stream.width(10);
            stream.precision(3);
            stream << node->Get_nDC(b);
        } else {
            stream.width(10);
            stream.precision(3);
            stream << 0.0;
        }
    }
    std::cout << stream.str() << std::endl;
}

void cGems::showSystem(Case& mel) {

    std::stringstream stream;

    std::cout << " TK : " << mel.TK << "   time : " << mel.Time << "   step : " << mel.i << std::endl;

    stream.str("");
    stream.clear();
    stream.fill(' ');    

    stream << "--------------Compounds----------------------" << std::endl;
    stream.width(12);
    stream << "DC";
    stream.width(8);
    stream << "index";
    stream.width(10);
    stream << "Moles" << std::endl;

    for (int i = 0; i < getDCompsCount(); i++) {
        long int b = node->DC_xCH_to_xDB(i);
        if (node->Get_nDC(b) > 1.0e-10) {
            stream.width(12);
            stream << getDCName(i);
            stream.width(8);
            stream.precision(3);
            stream << i;
            stream.width(10);
            stream.precision(3);
            stream << node->Get_nDC(b) << "\t";
            stream << node->Get_cDC(b) << std::endl;
        }
    }

    stream << "--------------Phases------------------------" << std::endl;
    stream.width(12);
    stream << "DC";
    stream.width(8);
    stream << "index";
    stream.width(10);
    stream << "Moles" << std::endl;

    for (int i = 0; i < getPhasesCount(); i++) {
        long int b = node->Ph_xCH_to_xDB(i);
        if (node->Ph_Mass(b) > 1.0e-20) {
            stream.width(12);
            stream << getPhaseName(i);
            stream.width(8);
            stream.precision(3);
            stream << i;
            stream.width(10);
            stream.precision(3);
            stream << node->Ph_Mass(b) << std::endl;
        }
    }

    std::cout << stream.str() << std::endl;
}



void cGems::readSetup(const std::string& file) {
    files.setup = file;

    using std::getline;
    ValueUnit ValUn;

    std::fstream setupFile;
    std::string line;
    vector<std::string> res;
    vector<string>::iterator begin, end;

    setup.clear();

    setupFile.open(file, std::ios::in);

    if (!setupFile.is_open()) {
        std::cout << "The setup data file [ " << file <<  " ] is not fould " << std::endl;
        return;
    }

    while (!setupFile.eof()) {
        std::getline(setupFile, line);

        res   = split(line);
        begin = res.begin();
        end   = res.end();

        if (line[0]!= '#' && line[0]!= '/' && !line.empty()) {

            // Read the data for the fit
            if ((begin != res.end()) && (*begin == "input" || *begin == "Input" || *begin == "INPUT")) {
                ++begin;
                while (begin != res.end()) {
                    setup.inputCols.push_back(*begin);
                    ++begin;
                }
                continue;
            }

            if ((begin != res.end()) && (*begin == "output" || *begin == "Output" || *begin == "OUTPUT"))
            {
                ++begin;
                while (begin != res.end()) {
                    setup.outputCols.push_back(*begin);
                    ++begin;
                }
                continue;
            }

            // Reads and initializes the gems project and the gems node
            if (*begin == "gems" || *begin == "Gems" || *begin == "GEMS") {
                ++begin;
                files.gems = *begin;
                initiateGems(files.gems);
                continue;
            }

            // add extra elements/compounds to the system
            if (*begin == "add" || *begin == "Add" || *begin == "ADD") {
                ++begin;
                std::string name = *begin;
                ++begin;
                data.add[name] = ToDouble(*begin);
                continue;
            }

            // Scaling
            if (*begin == "scale" || *begin == "Scale" || *begin == "SCALE") {
                ++begin;
                data.sc = ToDouble(*begin);
                continue;
            }

            // Maping compounds between MELCOR and GEMS
            if (*begin == "map" || *begin == "Map" || *begin == "MAP") {
                ++begin;
                std::string melName = *begin;
                ++begin;
                std::string gemName = *begin;
                data.nameMap[melName] = gemName;
                data.nameMap[gemName] = melName;
                continue;
            }


            // Printing the vapour pressure table
            if (*begin == "vcoeff" || *begin == "vtable") {
                ++begin;
                DataInfo::diPair pair;
                pair.str = *begin;
                ++begin;
                pair.i   = ToInt(*begin);
                data.pressureMap.push_back(pair);
                continue;
            }

            // Read the melcor file for the processing
            if (*begin == "meloutput" || *begin == "MELOUTPUT") {
                ++begin;
                files.melout = *begin;
                continue;
            }

            // Read the melcor file for the processing
            if (*begin == "meltable" || *begin == "melouttable") {
                ++begin;
                files.meltable = *begin;
                ++begin;
                while (begin != res.end()) {
                    setup.inputCols.push_back(*begin);
                    ++begin;
                }
                continue;
            }

            // Read the melcor file for the processing
            if (*begin == "melinput") {
                ++begin;
                files.melinp = *begin;
                continue;
            }

            // Read the time step for the melcor run
            if (*begin == "timestep") {
                ++begin;
                control.dTime = ToDouble(*begin);
                continue;
            }

            // Read the time step for the melcor run
            if (*begin == "time") {
                ++begin;
                control.Time = ToDouble(*begin);
                continue;
            }

            // Read the melcor file for the processing
            if (*begin == "read") {
                ++begin;
                // Read the melcor file for the processing
                if (*begin == "meloutput" || *begin == "MELOUTPUT") {
                    ++begin;
                    files.melout = *begin;
                    readMelcor();
                }
                continue;
            }

            // Read the melcor file for the processing
            if (*begin == "meloutput" || *begin == "MELOUTPUT") {
                ++begin;
                files.melout = *begin;
                readMelcor();
                continue;
            }

            // Read the melcor file for the processing
            if (*begin == "melcor" || *begin == "Melcor" || *begin == "MELCOR") {
                ++begin;
                files.melcor = *begin;
                continue;
            }

            // get the restart file name
            if (*begin == "restart" || *begin == "Restart" || *begin == "RESTART") {
                ++begin;
                files.restart = *begin;
                continue;
            }

            // Read the melcor file for the processing
            if (*begin == "melgen" || *begin == "Melgen" || *begin == "MELGEN") {
                ++begin;
                files.melgen = *begin;
                continue;
            }

            if (*begin == "ptf" || *begin == "Ptf" || *begin == "PTF") {
                ++begin;
                files.ptf = *begin;
                continue;
            }

            if (*begin == "population" || *begin == "Population" || *begin == "POPULATION") {
                ++begin;
                setup.population = ToInt(*begin);
                continue;
            }



            // RUN
            if (*begin == "run" || *begin == "RUN" || *begin == "Run") {
                std::string res = "";
                ++begin;

                if (*begin == "melcor" || *begin == "Melcor" || *begin == "MELCOR") {

                    std::string args = " ";

                    ++begin;
                    while (begin != end) {
                        args += *begin + " ";
                        ++begin;
                    }

                    res = files.melcor + args;
                    runMelcor(res);
                    continue;
                }

                if (*begin == "melgen" || *begin == "Melgen" || *begin == "MELGEN") {

                    std::string args = " ";

                    ++begin;
                    while (begin != end) {
                        args += *begin + " ";
                        ++begin;
                    }

                    res = files.melgen + args;
                    runMelgen(res);
                    continue;
                }

                if (*begin == "cgems" || *begin == "CGems" || *begin == "CGEMS") {
                    runCGEMS();
                    continue;
                }
            }


            // Read and set temperature of the calculations
            if (*begin == "T" || *begin == "temperature" || *begin == "Temperature") {
                ++begin;
                ValUn = getValue(*begin);
                if (ValUn.unut == "K")
                    setTemperatureK(ValUn.value);
                if (ValUn.unut == "C") {
                    setTemperatureC(ValUn.value-273.15);
                }
                continue;
            }

            if (*begin == "diffusion" || *begin == "Diffusion" || *begin == "DIFFUSION") {
                ++begin;
                std::string name = *begin;
                ++begin;
                data.diff[name].k0 = ToDouble(*begin);
                ++begin;
                data.diff[name].Q = ToDouble(*begin);
                continue;
            }

            // Read and set pressure of the calculations
            if (*begin == "P" || *begin == "pressure" || *begin == "Pressure") {
                ++begin;
                ValUn = getValue(*begin);
                if (ValUn.unut == "Bar" || ValUn.unut == "bar")
                    setPressureBar(ValUn.value);
                if (ValUn.unut == "Pa")
                    setPressurePa(ValUn.value);
                continue;
            }

        }
        //std::getline(setupFile, line);
    }

    setupFile.close();
    std::cout << "Done loading the setup file ... " << std::endl;
}

/*
void cGems::runPhysics(MelcorCase& mel) {
    if (mel.i > 0) {
        mel.core.elems    = melcor.cases[mel.i-1].core.elems;
        mel.surface.elems = melcor.cases[mel.i-1].surface.elems;
        mel.core.mol    = melcor.cases[mel.i-1].core.mol;
        mel.surface.mol = melcor.cases[mel.i-1].surface.mol;
    }

    diffuse(mel, mel.core, mel.surface);
    diffuse(mel, mel.gas,  mel.surface);
}
*/

void cGems::runMelcor(std::string& file) {
    std::string mel = "./melcor";
    std::system((mel + " " + file).c_str());
}

void cGems::runMelgen(std::string& file) {
    std::string mel = "./melgen";
    std::system((mel + " " + file).c_str());
}

void cGems::runCGEMS() {
    std::string res;
    readMelcor();
    std::cout << "Done ------ readMelcor" << std::endl;
    double time = 0;
    
    while (time < control.Time) {
        
        Case& cas = cases.back();
        std::cout << "Done ------ MelcorCase" << std::endl;
        
        //runPhysics(mel);
        std::cout << "Done ------ runPhysics" << std::endl;
        try {
            runAnalysis(cas);        // Run the gems equlibrium calculation
            std::cout << "Done ------ runAnalysis" << std::endl;
        } catch (std::runtime_error &err) { }
        
        getPressure(cas);
        std::cout << "Done ------ getPressure" << std::endl;

        cas.Time += control.dTime;
        time     += control.dTime;

        writeMelcorInput(cas);
        std::cout << "Done ------ writeMelcorInput" << std::endl;

        res = files.melinp + " " + "ow=e";
        runMelcor(res);
        std::cout << "Done ------ runMelcor" << std::endl;

        readMelcor();
        std::cout << "Done ------ readMelcor" << std::endl;
    }
}

void cGems::runGems(Case& mel) {
    // re-calculating equilibrium by calling GEMS3K, getting the status back
    dBR->NodeStatusCH = NEED_GEM_AIA;
    long int NodeStatusCH = node->GEM_run( false );
    if (NodeStatusCH == 2) mel.calc = true;
    else mel.calc = false;
    std::cout << "   NodeStatusCH  " << data.GemsStatus[NodeStatusCH] << std::endl;
}

void cGems::runAnalysis(Case& mel) {
    std::map<std::string, double>::const_iterator it;

    //setCase(mel);
    setCase(mel);

    runGems(mel);

    Compound compound;
    std::string name, phase, comp;
    std::map<std::string, int> decomp;
    std::map<std::string, double> elcompos;
    double mol, total, act;

    mel.compounds.clear();
    for (unsigned int i = 0; i < getDCompsCount(); i++) {
        Compound compound;
        long int b = node->DC_xCH_to_xDB(i);
        name = getDCName(i);
        mol  = node->Get_nDC(b);
        act  = node->Get_aDC(b, false);
        phase= getPhase(name);
        comp = cleanName(name);

        elcompos = getElemCompos(comp, mol);

        total = 0.0;
        for (it = elcompos.begin(); it != elcompos.end(); it++) {
            compound.elems[it->first]  += it->second;
            total += it->second * data.molMass[it->first] * 0.001;
        }
        compound.mol   = mol;
        compound.kg    = total;
        compound.name  = comp;
        compound.phase = phase;
        compound.prop.act   = act;
        compound.prop.P = pow(10.0, act) * 750.062;
        mel.compounds.push_back(compound);
    }
}
