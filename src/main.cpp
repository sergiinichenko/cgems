//--------------------------------------------------------------------
// $Id: main.cpp 686 2012-06-01 14:10:22Z kulik $
//
//
//-------------------------------------------------------------------

#include "cgems.h"

int main(int argc, char *argv[]) {


    cGems cgems;

    if (argc > 1) {
        cgems.files.setup = argv[1];
    }

    //cgems.files.setup = "setup.inp";
    cgems.readSetup(cgems.files.setup);

    /*
    std::string comp;
    std::map<std::string, double> sample;
    std::map<std::string, double> elcomp;
    std::map<std::string, double>::const_iterator it;
    std::map<std::string, double>::const_iterator el;
    std::map<std::string, int>::const_iterator form;

    sample["CsI"] = 2.0;
    sample["CsI(g)"] = 2.0;
    sample["I(g)"] = 2.0;
    sample["I2(g)"] = 2.0;
    sample["H2O"] = 2.0;
    sample[""] = 2.0;
    sample["CsI"] = 2.0;
    sample["CsI"] = 2.0;
    sample["CsI"] = 2.0;

    for (it = sample.begin(); it != sample.end(); it++) {
        comp = cgems.cleanName(it->first);
        std::map<std::string, int> res = cgems.parser.readFormula(comp);
        std::cout << it->first << " : " << it->second << " - " << comp << std::endl;

        for (form = res.begin(); form != res.end(); form++) {
            std::cout << "   " << form->first << " : " <<
form->second << std::endl;
        }
        elcomp = cgems.getElemCompos(comp, it->second);
        for (el = elcomp.begin(); el != elcomp.end(); el++) {
            std::cout << "       " << el->first << " : " << el->second << std::endl;
        }

    }
    */

    //cgems.rePopulateData();

    /*
    for (MelcorCase& mel : cgems.melcor.cases) {
        cgems.runPhysics(mel);
        cgems.runAnalysis(mel);        // Run the gems equlibrium calculation
        cgems.getPressure(mel);
        cgems.removeGas(mel);
        cgems.showSystem(mel);     // Print results to the screen
        cgems.writeMelcorInput(mel);
    }
    cgems.saveDF();
    */


    /*
    cgems.setTemperatureC(1700.0);
    cgems.pauseMELCOR();
    cgems.readPTF(ptffile);
    */
    return 0;
}
