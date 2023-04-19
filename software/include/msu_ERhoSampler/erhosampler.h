#ifndef __MSU_ERHO_SAMPLER_H__
#define __MSU_ERHO_SAMPLER_H__
#include <Eigen/Dense>
#include "msu_eos/resonances.h"
#include "msu_sampler/part.h"
#include "msu_sampler/classdefs.h"
#include "msu_sampler/master.h"
#include "msu_sampler/sampler.h"
#include "msu_eos/eos.h"
#include "msu_commonutils/log.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/decay_nbody.h"
using namespace std;

namespace NMSU_ERrhoSampler{
	void FillOutHyperBjorken(Chyper *hyper,double T,double tau,double R,double deleta,double rhoB,double rhoII);
	void DecayParts(Crandy *randyset,CpartList *partlist);
	void IncrementEQtest(CpartList *partlist,Eigen::VectorXd &EQtot,Eigen::VectorXd &EQTarget);
	void Chi4Test(CpartList *partlist,Eigen::MatrixXd &chitest);
	void GetDecayCorrs(CparameterMap *parmap,Crandy *randy,CpartList *motherpartlist);
}

class CcorrVsEta{
public:
	double DETA;
	virtual void GetCorrVsEta(double eta,Eigen::MatrixXd &corr){
		//This is a dummy
	}
};
class CcorrVsEtaScott : public CcorrVsEta{
public:
	void GetCorrVsEta(double eta,Eigen::MatrixXd &corr);
};

class CcorrVsEtaOleh : public CcorrVsEta{
public:
	void GetCorrVsEta(double eta,Eigen::MatrixXd &corr);
};

class CcorrVsY{
public:
	CparameterMap *parmap;
	double DY;
	int NY;
	Crandy *randy;
	CcorrVsY(CparameterMap *parmap);
	vector<Eigen::MatrixXd> corr;
	double denom;
	void Increment(CpartList *partlista,CpartList *partlistb,CcorrVsEta *corrvseta);
	double GetPairWeight(Cpart *part1,Cpart *part2,Eigen::MatrixXd &CorrMatrix);
	void WriteResults();
	
};

class CdecayCorrVsY{
public:
	CparameterMap *parmap;
	double DY;
	int NY;
	CdecayCorrVsY(CparameterMap *parmap);
	vector<Eigen::MatrixXd> corr;
	void Increment(CpartList *partlist);
	void WriteResults();
};




#endif
