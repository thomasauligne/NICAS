//---------------------------------------------------------------------
/// Purpose: random number generator class implementation
/// Author : Benjamin Menetrier
/// Licensing: this code is distributed under the CeCILL-C license
/// Copyright © 2017 METEO-FRANCE
// ----------------------------------------------------------------------
#include "type_randgen.h"
#include "type_randgen.hpp"
#include "external/Cover_Tree.h"
#include "external/Cover_Tree_Point.h"
#include <ostream>
#include <iomanip>
#include <cmath>
#include <random>

using namespace std;

// Constructor
randGen::randGen(int default_seed) {
    // Initialize random number generator
    if (default_seed==0) {
        std::random_device rd;
        gen = new std::mt19937(rd());
    }
    else {
        gen = new std::mt19937;
    }
}

// Destructor
randGen::~randGen(){}

// Random integer generator
void randGen::rand_integer(int binf, int bsup, int *ir) const{
    // Initialize uniform distribution
    std::uniform_int_distribution<int> dis(binf,bsup);

    // Generate random integer
    *ir=dis(*gen);
    return;
}

// Sampling initialization
void randGen::initialize_sampling(int ddis, int n, double lon[], double lat[], int mask[], double L[], int ntry, int nrep, int ns, int nfor, int ifor[], int ihor[]) const{
    // Initialize uniform distribution
    std::uniform_int_distribution<int> dis(0,n-1);

    // Declaration
    int progint;
    int ir;

    // Copy mask (updated then)
    int mask_copy[n];
    for(int i=0;i<n;i++) {
        mask_copy[i]=mask[i];
    }

    // Initialize tree
    std::cout << std::string(10,' ') << "Compute sampling: ";
    CoverTree<CoverTreePoint> cTree(1.0e10);
    int is=0;

    // Include forced points
    for(int i=0;i<nfor;i++) {
        ihor[is]=ifor[i];
        cTree.insert(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]));
        mask_copy[ifor[i]-1]=0;
        is++;
    }

    // Fill the tree
    progint=ddis;
    std::cout << " 0% " << std::flush;
    while(is<ns) {
        // Find new point
        double distmax=0.0;
        int irmax=-1;
        int itry=0;
        while(itry<ntry) {
            // Generate random number
            ir=dis(*gen);
            if(mask_copy[ir]==1) {
                if(is>0) {
                    // Find nearest neighbor
                    vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(-999,lon[ir],lat[ir]),1));

                    // Check distance
                    double dist=neighbors[0].getDist()/sqrt(0.5*(pow(L[neighbors[0].getIndex()],2.0)+pow(L[ir],2.0)));
                    if(dist>distmax) {
                        distmax=dist;
                        irmax=ir;
                    }
                }
                else {
                    irmax=ir;
                }
            }
            itry++;
        }

        // Insert point
        if(irmax>0) {
            ihor[is]=irmax+1;
            cTree.insert(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]));
            mask_copy[irmax]=0;
            is++;
        }

        // Print progression
        double prog=100.0*(float)(is)/(float)ns;
        if((int)(prog)>progint) {
           std::cout << " " << (int)(progint) << "% " << std::flush;
           progint=progint+ddis;
        }
    }
    std::cout << "100%" << std::endl;

    // Get minimum distance
    double distmininit=HUGE_VAL;
    for(int is=nfor;is<ns;is++) {
        vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]),2));
        double dist=neighbors[1].getDist()/sqrt(0.5*(pow(L[neighbors[1].getIndex()],2.0)+pow(L[ihor[is]-1],2.0)));
        if(dist<distmininit) {
            distmininit=dist;
        }
    }

    if (nrep>0) {
        // Improve sampling with replacements
        std::cout << std::string(10,' ') << "Improve sampling: ";
        int irep=0;
        progint=ddis;
        std::cout << " 0% " << std::flush;
        while(irep<nrep) {
            // Get minimum distance
            double distmin=HUGE_VAL;
            int ismin=0;
            for(int is=nfor;is<ns;is++) {
                vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]),2));
                double dist=neighbors[1].getDist()/sqrt(0.5*(pow(L[neighbors[1].getIndex()],2.0)+pow(L[ihor[is]-1],2.0)));
                if(dist<distmin) {
                    distmin=dist;
                    ismin=is;
                }
            }

            // Remove point
            cTree.remove(CoverTreePoint(ismin,lon[ihor[ismin]-1],lat[ihor[ismin]-1]));

            // Find new point
            double distmax=0.0;
            int irmax=-1;
            int itry=0;
            while(itry<ntry) {
                // Generate random number
                ir=dis(*gen);
                if(mask_copy[ir]==1) {
                    // Find nearest neighbor
                    vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(-999,lon[ir],lat[ir]),1));

                    // Check distance
                    double dist=neighbors[0].getDist()/sqrt(0.5*(pow(L[neighbors[0].getIndex()],2.0)+pow(L[ir],2.0)));
                    if((dist>distmax) && (dist>distmininit)) {
                        distmax=dist;
                        irmax=ir;
                    }
                }
                itry++;
            }

            // Insert point
            if(irmax>0) {
                // Insert new point
                ihor[ismin]=irmax+1;
                cTree.insert(CoverTreePoint(ismin,lon[ihor[ismin]-1],lat[ihor[ismin]-1]));
                mask_copy[irmax]=0;
            }
            else {
                // Re-insert old point
                cTree.insert(CoverTreePoint(ismin,lon[ihor[ismin]-1],lat[ihor[ismin]-1]));
            }
            irep++;

            // Print progression
            double prog=100.0*(float)(irep)/(float)nrep;
            if((int)(prog)>progint) {
               std::cout << " " << (int)(progint) << "% " << std::flush;
               progint=progint+ddis;
            }
        }
        std::cout << "100%" << std::endl;

        // Find final minimum distance
        double distmin=HUGE_VAL;
        for(int is=nfor;is<ns;is++) {
            vector<CoverTreePoint> neighbors(cTree.kNearestNeighbors(CoverTreePoint(is,lon[ihor[is]-1],lat[ihor[is]-1]),2));
            double dist=neighbors[1].getDist()/sqrt(0.5*(pow(L[neighbors[1].getIndex()],2.0)+pow(L[ihor[is]-1],2.0)));
            if(dist<distmin) {
                distmin=dist;
            }
        }
        std::cout << std::string(10,' ') << "Minimum normalized separation increased from " << std::fixed << std::setprecision(5) << distmininit << " to " << distmin << std::endl;
    }
    else {
        std::cout << std::string(10,' ') << "Minimum normalized separation " << std::fixed << std::setprecision(5) << distmininit << std::endl;
    }
    return;
}

