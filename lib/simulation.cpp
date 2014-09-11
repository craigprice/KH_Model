#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
const static double PI = 3.14159265358979323846;
#include "ClassicalSpin3D.h"
#include "MonteCarlo.h"
#include "ClassicalSpin3D.cpp"
#include "MonteCarlo.cpp"
using namespace std;
int main(int argc, char*argv[]){
    if(argc != 27 && argc != 4){
        std::cerr << "Need inputs" << std::endl;
        return 1;
    }
    //Set the random seed.
    srand48( (unsigned long) time(0) );
    MonteCarlo* MC;// = new MonteCarlo();
    //delete MC;
    if(argc == 27){
        bool    isDomainWalls = (bool) (atoi(argv[2]));
        int     directionOfDomain =     atoi(argv[3]);
        bool    isRotHam    = (bool) (  atoi(argv[4]));

        bool    isAlphaHam  = (bool) (  atoi(argv[5]));
        double  alpha       =           atof(argv[6]);

        bool    isJ23Term   = (bool) (  atoi(argv[7]));
        double  j2          =           atof(argv[8]);
        double  j3          =           atof(argv[9]);
        double  k2          =           atof(argv[10]);

        bool    isCubeHam   = (bool) (  atoi(argv[11]));
        double  cubicD      =           atof(argv[12]);
        
        // B field Not necessarily normalized
        bool    isBField    = (bool) (  atoi(argv[13]));
        double  bField_x    =           atof(argv[14]);
        double  bField_y    =           atof(argv[15]);
        double  bField_z    =           atof(argv[16]);
        double  bFieldMag   =           atof(argv[17]);
        
        bool    isPhiHam    = (bool) (  atoi(argv[18]));
        double  phi         =           atof(argv[19]);
        
        int     cellsA      =           atoi(argv[20]);
        int     cellsB      =           atoi(argv[21]);
        int     cellsC      =           atoi(argv[22]);
        
        double  KbT         =           atof(argv[23]);

        double  estimatedTc =           atof(argv[24]);
        int     numSweepsToPerformTotal = atoi(argv[25]);
        std::string path    =           argv[26];
        MC = new MonteCarlo(isDomainWalls,
                            directionOfDomain,
                            isRotHam,
                            isAlphaHam,
                            alpha,
                            isJ23Term,
                            j2,
                            j3,
                            k2,
                            isCubeHam,
                            cubicD,
                            isBField,
                            bField_x,
                            bField_y,
                            bField_z,
                            bFieldMag,
                            isPhiHam,
                            phi,
                            cellsA,
                            cellsB,
                            cellsC,
                            KbT,
                            estimatedTc,
                            numSweepsToPerformTotal,
                            path);
    }else if(argc == 4){
        bool reThermalize = (bool) (atoi(argv[2]));
        std::string filename = argv[3];
        MC = new MonteCarlo(reThermalize, filename);
    }else{
        cerr<<"inputs"<<endl;
        exit(3161990);
    }
    //int numSweeps             =  500 * 1000;
    double durationOfSingleRun = 0.35;
    int minSweepsToThermalize = 100 * 1000;
    int numSweepsBetwnConfigs = 6;
    bool outOfTime = false;
    
    //
    //Runs the lattice without taking data at these temperatures to get the
    //spins pointed in approximately the correct direction.
    
    if(MC -> MD.numSweepsUsedToThermalize < (int) (minSweepsToThermalize / 2)){
        outOfTime = MC -> thermalize(3.0 * MC -> MCP.estimatedTc,
                                     (int) (minSweepsToThermalize / 2),
                                     durationOfSingleRun);
    }
    if((!outOfTime)&&
       (MC -> MD.numSweepsUsedToThermalize < minSweepsToThermalize)){
        outOfTime = MC -> thermalize(MC -> MCP.KbT,
                                     minSweepsToThermalize,
                                     durationOfSingleRun);
    }
    
    while (!(MC -> MD.isThermal) && !outOfTime) {
        outOfTime = MC -> isThermalized(durationOfSingleRun);
    }
    
    if(!outOfTime && MC -> MD.isThermal){
        //To run a snapshot, change numSweeps, sweep().
        //MC -> sweepSnapshot(numSweeps, 6, 9);//This hasn't been edited - check function
        
        MC -> sweep(numSweepsBetwnConfigs,
                    durationOfSingleRun);
    }
    
    MC->finalPrint();
    std::cout << "done!" << std::endl;
    return 0;
}
