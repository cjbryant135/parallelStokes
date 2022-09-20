/*
 * Driver for solving the Stokes equations on a polar grid. Sets up a Stokes momentum problem
 * for a given pressure/body forces and solves via SOR. 
 *	
 *		Solution (with respect to polar coordinates):
 *			u = r^2 cos(th)
 *			v = -3 r^2 sin(th)
 *			p = r sin(th)
 *
 *		Body force:
 *			f_r = 8 cos(th)-sin(th)
 *			f_th = -8 sin(th)-cos(th) 
 *		
 *		Inputs:
 *			nr - number of grid points in r-direction 
 *			nt - number of grid points in t-direction 
 *			outDir - location to write data to e.g. "./" Assumes this directory already exists. 
 *	
 *		Outputs:
 *			-Writes number of SOR iterations needed to the console. 
 *			-Flow data written to "<outDir>u.out", "<outDir>v.out", "<outDir>p.out"
 */

//includes
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<string>

#include "polarStokesGrid.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	
	//user inputs 
	int nr = atoi(argv[1]); //number of grid points in r direction
	int nt = atoi(argv[2]); //number of grid points in theta direction 
	string outDir(argv[3]); //directory to write data to (assumes this directory exists)

	//set the domain
	double rP = 0.1;
	double rMax = 0.4;
	
	//create the grid 
	polarStokesGrid annulus(nr, nt, rP, rMax);

	//set numerical parameters
	annulus.setSORTol(1.e-8); //tolerance for the SOR solver
	annulus.setOmega(1.8); //relaxation parameter
	
	//set boundary conditions
	int i,j;
	double rG, tG;
	for(j=0; j<nt-1; ++j) {
		//r-vel BC (note: the SOR solver does not touch u values at i=0 or i=nr-1)
		tG = annulus.tC(j);
		annulus.setU(nr-1,j,rMax*rMax*cos(tG));
		annulus.setU(0,j,rP*rP*cos(tG));

		//t-vel BC
		tG = annulus.tF(j);
		annulus.setVRMax(j,-3.*rMax*rMax*sin(tG));
		annulus.setVRMin(j,-3.*rP*rP*sin(tG));
	}
	
	//set body forces
	for(i=0; i<nr; ++i) {
		for(j=0; j<nt-1; ++j) {
			tG = annulus.tC(j);
			annulus.setFR(i,j,8.*cos(tG)-sin(tG));
			if(i<nr-1) {
				tG = annulus.tF(j);
				annulus.setFT(i,j,-8.*sin(tG)-cos(tG));
			}
		}
	}

	//set the correct pressure
	for(i=0; i<nr-1; ++i) {
		for(j=0; j<nt-1; ++j) {
			rG = annulus.rC(i);
			tG = annulus.tC(j);
			annulus.setP(i,j,rG*sin(tG));	
		}
	}


	//solve Stokes momentum equations
	int count = annulus.solveMomSOR();
	if(count < 1e5)
		cout << "Stokes momentum equation solved in " << count << " iterations\n";
	else {
		cout << "ERROR: Solver did not converge in 100000 iterations. Exiting...\n";
		exit(1);
	}
	//dump data and clean up
	string uOut = outDir + "u.out";
	string vOut = outDir + "v.out";
	string pOut = outDir + "p.out";

	annulus.dumpFlowData(uOut,vOut,pOut);

	//exit success 
	return 0;
}
