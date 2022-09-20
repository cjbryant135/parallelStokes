#include "polarStokesGrid.hpp"

//definitions for the polarStokesGrid class

//****one-liners****
//--setters--
void polarStokesGrid::setOmega(const double val) {omega=val;} //sets relaxation parameter 
void polarStokesGrid::setSORTol(const double val) {sortol=val;} //sets tolerance for SOR solver
void polarStokesGrid::setVRMin(const int j, const double val) {vRMin[j]=val;} //v-velocity data at r=rMin
void polarStokesGrid::setVRMax(const int j, const double val) {vRMax[j]=val;} //v-velociy data at r=rMax
void polarStokesGrid::setU(const int i, const int j, const double val) {flowData[uI(i,j)]=val;} //set u-velocity at any point
void polarStokesGrid::setV(const int i, const int j, const double val) {flowData[vI(i,j)]=val;} //set v-velocity at any point 
void polarStokesGrid::setP(const int i, const int j, const double val) {flowData[pI(i,j)]=val;} //set pressure data at any pount 
void polarStokesGrid::setFR(const int i, const int j, const double val) {bodyRData[uI(i,j)]=val;} //set r body force
void polarStokesGrid::setFT(const int i, const int j, const double val) {bodyTData[vI(i,j)-vStart]=val;} //set theta body force 
void polarStokesGrid::turnOnVerbose() {verbose=true;}; //writes residual information during the solve 

//--data access--
//grid sizes 
int polarStokesGrid::Nr() const {return mNr;} 
int polarStokesGrid::Nt() const {return mNt;}
//flow data and body forces 
double polarStokesGrid::u(const int i, const int j) const {return flowData[uI(i,j)];}
double polarStokesGrid::v(const int i, const int j) const {return flowData[vI(i,j)];}
double polarStokesGrid::p(const int i, const int j) const {return flowData[pI(i,j)];}
double polarStokesGrid::fr(const int i, const int j) const {return bodyRData[uI(i,j)];}
double polarStokesGrid::ft(const int i, const int j) const {return bodyTData[vI(i,j)-vStart];}
//grid point physical locations
double polarStokesGrid::rC(const int i) const {return mRMin+(i+0.5)*dr;}
double polarStokesGrid::rF(const int i) const {return mRMin+i*dr;}
double polarStokesGrid::tC(const int j) const {return (j+0.5)*dt;}
double polarStokesGrid::tF(const int j) const {return j*dt;}
double polarStokesGrid::xULoc(const int i, const int j) const {return rF(i)*cos(tC(j));}
double polarStokesGrid::yULoc(const int i, const int j) const {return rF(i)*sin(tC(j));}
double polarStokesGrid::xVLoc(const int i, const int j) const {return rC(i)*cos(tF(j));}
double polarStokesGrid::yVLoc(const int i, const int j) const {return rC(i)*sin(tF(j));}
double polarStokesGrid::xPLoc(const int i, const int j) const {return rC(i)*cos(tC(j));}
double polarStokesGrid::yPLoc(const int i, const int j) const {return rC(i)*sin(tC(j));}
//grid point indices in flowdata array 
int polarStokesGrid::uI(const int i, const int j) const {return i*(mNt-1)+j;}
int polarStokesGrid::vI(const int i, const int j) const {return vStart+i*(mNt-1)+j;}
int polarStokesGrid::pI(const int i, const int j) const {return pStart+i*(mNt-1)+j;}


//****constructor****
polarStokesGrid::polarStokesGrid(const int nr, const int nt, const double rMin, const double rMax) {
	//set parameters
	mNr = nr;
	mNt = nt;
	mRMin = rMin;
	mRMax = rMax;

	//calculate dr/dt
	dr = (rMax-rMin)/(nr-1);
	dt = 2.*M_PI/(nt-1);

	//set default numerical parameter values 
	sortol = 1.e-6;
	omega = 1.8;
	
	//allocate data arrays
	flowData = new double[mNr*(mNt-1)+(mNr-1)*(mNt-1)+(mNr-1)*(mNt-1)]();
	bodyRData = new double[mNr*(mNt-1)]();
	bodyTData = new double[(mNr-1)*(mNt-1)]();
	vRMax = new double[mNt-1]();
	vRMin = new double[mNt-1]();

	//set starting indices for the flow data 
	uStart = 0;
	vStart = mNr*(mNt-1);
	pStart = vStart+(mNr-1)*(mNt-1);

	verbose = false; 
}

//solver routines 
void polarStokesGrid::makeResidualsHuge() {
	maxURes = 1.e10;
	maxVRes = 1.e10;
	maxPRes = 1.e10;
}

void polarStokesGrid::zeroOutResiduals() {
	maxURes = 0.;
	maxVRes = 0.;
	maxPRes = 0.;
}

bool polarStokesGrid::isConverged() {
	return maxURes < sortol && maxVRes < sortol && maxPRes < sortol;
}

int polarStokesGrid::solveMomSOR() {
	if(verbose) cout << "==============================================================\n";
	//storage 
	int count = 0;
	int i,j;
	int jP, jM;
	double here, west, east, nort, sout;
	double tempRes;
	//compute prefactors 
	double drDt = dr/dt;
	double dtDr = dt/dr;
	makeResidualsHuge(); //force us into the loop
	while(!isConverged() && count < MAXITS) {
		count++;
		if(verbose && count%100==0) {
			cout << "On iteration " << count << " uRes = " << maxURes << ", vRes = " << maxVRes << endl;
		}
		zeroOutResiduals();
		//update r velocity
		for(i=1; i<mNr-1; ++i) { //note: skip grid points at particle surface and r=rmax
			for(j=0; j<mNt-1; ++j) {
				jM = j-1;
				jP = j+1;
				//periodicity 
				if(jM < 0) jM = mNt-2;
				if(jP > mNt-2) jP = 0;
				here = u(i,j);
				nort = u(i,jP);
				sout = u(i,jM);
				east = u(i+1,j);
				west = u(i-1,j);

				tempRes = dtDr * rF(i) * (east - 2.*here + west) + 
							drDt * (nort - 2.*here + sout)/rF(i) + 
							dt * (east - west)/2. - 
							dr * ((v(i-1,jP)+v(i,jP))-(v(i-1,j)+v(i,j)))/rF(i) - 
							dr * dt * here/rF(i) - 
							dt * rF(i) * (p(i,j)-p(i-1,j)) - 
							dr * dt * rF(i) * bodyRData[uI(i,j)-uStart];
				//update
				flowData[uI(i,j)] += omega*tempRes/(dr*dt/rF(i) + 2.*(rF(i)*dtDr+drDt/rF(i)));
				if( fabs(tempRes) > maxURes ) maxURes = fabs(tempRes);
			}
		}

		//update t velocity
		for(i=0; i<mNr-1; ++i) { 
			for(j=0; j<mNt-1; ++j) {
				jM = j-1;
				jP = j+1;
				//periodicity 
				if(jM < 0) jM = mNt-2;
				if(jP > mNt-2) jP = 0;
				here = v(i,j);
				nort = v(i,jP);
				sout = v(i,jM);
				//cases for inner/outer boundaries
				if(i==0) {
					west = 2.*vRMin[j] - here;
				} else {
					west = v(i-1,j);
				}
				if(i==mNr-2) {
					east = 2.*vRMax[j] - here;
				} else {
					east = v(i+1,j);
				}
				
				tempRes = dtDr * rC(i) * ( east - 2.*here + west ) + 
							drDt * (nort - 2.*here + sout)/rC(i) + 
							dt * (east-west)/2. + 
							dr * ((u(i,j)+u(i+1,j))-(u(i,jM)+u(i+1,jM)))/rC(i) - 
							dr * dt * here/rC(i) - 
							dr * (p(i,j)-p(i,jM)) - 
							dr * dt * rC(i) * bodyTData[vI(i,j)-vStart];
				//update	
				flowData[vI(i,j)] += omega*tempRes/(dr*dt/rC(i) + 2*(rC(i)*dtDr + drDt/rC(i)));
				if( fabs(tempRes) > maxVRes ) maxVRes = fabs(tempRes);

			}
		}

	}
	if(verbose) cout << "==============================================================\n";
	return count;	

}

void polarStokesGrid::dumpFlowData(string uOut, string vOut, string pOut) {
	FILE* uFile;
	FILE* vFile;
	FILE* pFile;
	uFile = fopen(uOut.c_str(), "w");
	vFile = fopen(vOut.c_str(), "w");
	pFile = fopen(pOut.c_str(), "w");
	
	int urLen, utLen, vrLen, vtLen, prLen, ptLen;
	urLen = mNr;
	utLen = mNt-1;
	vrLen = mNr-1;
	vtLen = mNt-1;
	prLen = mNr-1;
	ptLen = mNt-1;
	
	//write array sizes
	fwrite(&urLen, sizeof(int), 1, uFile);
	fwrite(&utLen, sizeof(int), 1, uFile);
	fwrite(&vrLen, sizeof(int), 1, vFile);
	fwrite(&vtLen, sizeof(int), 1, vFile);
	fwrite(&prLen, sizeof(int), 1, pFile);
	fwrite(&ptLen, sizeof(int), 1, pFile);
	
	//write data
	fwrite(&(flowData[uStart]), sizeof(double), urLen*utLen, uFile);
	fwrite(&(flowData[vStart]), sizeof(double), vrLen*vtLen, vFile);
	fwrite(&(flowData[pStart]), sizeof(double), prLen*ptLen, pFile);
	
	//create coordinate arrays
	double* uXCoords = new double[urLen*utLen];
	double* uYCoords = new double[urLen*utLen];
	double* vXCoords = new double[vrLen*vtLen];
	double* vYCoords = new double[vrLen*vtLen];
	double* pXCoords = new double[prLen*ptLen];
	double* pYCoords = new double[prLen*ptLen];
	
	int i,j;
	for(i=0; i<urLen; ++i) {
		for(j=0; j<utLen; ++j) {
			uXCoords[i*utLen+j] = xULoc(i,j);
			uYCoords[i*utLen+j] = yULoc(i,j);
		}
	}
	for(i=0; i<vrLen; ++i) {
		for(j=0; j<vtLen; ++j) {
			vXCoords[i*vtLen+j] = xVLoc(i,j);
			vYCoords[i*vtLen+j] = yVLoc(i,j);
		}
	}
	for(i=0; i<prLen; ++i) {
		for(j=0; j<ptLen; ++j) {
			pXCoords[i*ptLen+j] = xPLoc(i,j);
			pYCoords[i*ptLen+j] = yPLoc(i,j);
		}
	}
	
	//write coordinates
	fwrite(uXCoords, sizeof(double), urLen*utLen, uFile);
	fwrite(uYCoords, sizeof(double), urLen*utLen, uFile);
	fwrite(vXCoords, sizeof(double), vrLen*vtLen, vFile);
	fwrite(vYCoords, sizeof(double), vrLen*vtLen, vFile);
	fwrite(pXCoords, sizeof(double), prLen*ptLen, pFile);
	fwrite(pYCoords, sizeof(double), prLen*ptLen, pFile);
	//clean up
	delete[] uXCoords;
	delete[] uYCoords;
	delete[] vXCoords;
	delete[] vYCoords;
	delete[] pXCoords;
	delete[] pYCoords;


	fclose(uFile);
	fclose(vFile);
	fclose(pFile);
}

//destructor
polarStokesGrid::~polarStokesGrid(void) {
	if(flowData) delete[] flowData;
	if(bodyRData) delete[] bodyRData;
	if(bodyTData) delete[] bodyTData;
	if(vRMax) delete[] vRMax;
	if(vRMin) delete[] vRMin;
}
