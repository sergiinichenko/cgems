#ifndef READER_H
#define READER_H

#include "datatypes.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include "formulaparser.h"
#include "physics.h"
#include "eigen/Core"
#include "eigen/Dense"

class Reader {
public:
    DATABR* dBR;
    TNode* node;


    SetupStruct setup;
    DataInfo data;
    FilesStruct files;
    Controls control;
    std::vector<Case> cases;

    FormulaParser parser;

    bool fileExists(const std::string& filename);
    std::vector<std::string> split(const std::string& str);
    ValueUnit getValue(std::string& str);
    double ToDouble(const std::string&);
    int ToInt(const string& s);

    unsigned long getDCompsCount();
    unsigned long getICompsCount();
    unsigned long getPhasesCount();
    std::string getDCName(unsigned int);
    std::string getICName(unsigned int);
    std::string getPhaseName(int);
    std::string getPhase(const std::string&);
    std::string cleanName(const std::string&);
    void setAmmountOfDC(std::string name, double val);
    void setAmmountOfIC(std::string name, double val);
    void setTemperatureK(double val);
    void setTemperatureC(double val);
    void setPressureBar(double);
    void setPressurePa(double);
    void readMelcor();
    void readMelcorOutput(const std::string&);
    void readMelcorTable(const std::string&);
    void pauseMELCOR();
    void pauseMELCOR(const std::string&);
    void setCase(Case&);
    double getMolarMass(const std::string&);
    std::map<std::string, double> getElemCompos(const std::string&, const double&);
    std::map<std::string, double> getElemActivities(const std::string&, const double&);
    void getPressure(Case&);
    std::vector<double> CCPressureCoeffs(std::string&);
    double getTboil(std::string comp);
    double getTmelt(std::string comp);
    double elementChemPotential(std::string&);
    double dcActivity(std::string&);
    void getElements(Case&);
    void saveDF();
    bool findNextLine(std::fstream& file, std::string& find, std::string& line);
    bool findNextLinePart(std::fstream& file, std::string& find, std::string& line);
    bool isEqual(std::string& find, std::string& line);
    void diffuse();
    void diffuse(Case&);
    void diffuse(Case& from, Case& to);
    void rePopulateData();
    std::map<std::string, double> elements();
    void printMELCOR(Case&, std::string);
    void printGEMS(Case&, std::string);
    void addStrToVec(std::vector<std::string>&, std::string);
    void writeMelcorInput(Case&);
    std::vector<double> calcCoeffs(std::vector<std::vector<double>> x, std::vector<double> y);
    void removeGas(Case&);
    void clearFiles(std::string& f);
    int getIndex(std::vector<Compound>&, std::string&, std::string&);
    int getIndex(std::vector<Compound>&, std::string&);
    bool contains(std::string, std::string);

    Reader();
    //~Reader();
};

#endif // READER_H
