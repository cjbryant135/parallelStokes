/*
 * Simplified version of the polarStokesGrid class used in the overset code. Uses SOR to solve
 * the Stokes momentum equation for a given pressure field.
 *
 */

#ifndef POLARSTOKESGRID_H
#define POLARSTOKESGRID_H

//includes 
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<algorithm>

using namespace std;

class polarStokesGrid
{
	static const int MAXITS = 1e5; //max iterations for the SOR solve
	public: 
		//constructor 
		polarStokesGrid(const int nr, const int nt, const double rMin, const double rMax);
		
		//setters 
		void setVRMin(const int j, const double val);
		void setVRMax(const int j, const double val);
		void setU(const int i, const int j, const double val);
		void setV(const int i, const int j, const double val);
		void setP(const int i, const int j, const double val);
		void setFR(const int i, const int j, const double val);
		void setFT(const int i, const int j, const double val);
		void setOmega(const double val);
		void setSORTol(const double val);

		//data access 
		int Nr() const;
		int Nt() const;
		double u(const int i, const int j) const;
		double v(const int i, const int j) const;
		double p(const int i, const int j) const;
		double fr(const int i, const int j) const;
		double ft(const int i, const int j) const;

		
		//grid point locations/indices
		double rC(const int i) const;
		double tC(const int j) const;
		double rF(const int i) const;
		double tF(const int j) const;
		double xULoc(const int i, const int j) const;
		double yULoc(const int i, const int j) const;
		double xVLoc(const int i, const int j) const;
		double yVLoc(const int i, const int j) const;
		double xPLoc(const int i, const int j) const;
		double yPLoc(const int i, const int j) const;
		
		int uI(const int i, const int j) const;
		int vI(const int i, const int j) const;
		int pI(const int i, const int j) const;

		//solver routines 
		void makeResidualsHuge();
		void zeroOutResiduals();
		bool isConverged();

		int solveMomSOR();

		void dumpFlowData(string uOut, string vOut, string pOut);
		
		~polarStokesGrid(void);
	private:
		double* flowData;
		double* bodyRData;
		double* bodyTData;
		double* vRMin;
		double* vRMax;

		double mRMin, mRMax;
		double dr, dt;
		int mNr, mNt;

		double sortol;
		double omega;
		double maxURes, maxVRes, maxPRes;

		int uStart, vStart, pStart;
};


#endif
