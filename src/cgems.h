#ifndef CGEMS_H
#define CGEMS_H

#include "reader.h"
#include "node.h"

class cGems : public Reader {
public:
    std::string path;

    cGems();
    cGems(std::string file);
    void initiateGems(std::string file);
    void showSystem(Case&);
    void showSelected();
    void controlFile(std::string);
    void readSetup(const std::string&);

    void runPhysics();
    void runPhysics(Case&);
    void runMelcor(std::string&);
    void runMelgen(std::string&);
    void runCGEMS();
    void runGems(Case&);
    void runAnalysis(Case&);

};

#endif // CGEMS_H
