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
using namespace std;

class CcorrVsEta{
public:
	double DETA;
	void GetCorrVsEta(double eta,Eigen::MatrixXd &corr);
};

class CcorrVsY{
public:
	double DY;
	int NY;
	Crandy *randy;
	CcorrVsY(CparameterMap *parmap);
	vector<Eigen::MatrixXd> corr;
	vector<double> denom;
	void Increment(CpartList *partlista,CpartList *partlistb,CcorrVsEta *corrvseta);
	GetPairWeight(Cpart *part1,Cpart *part2,Eigen::MatrixXd &CorrMatrix);
	void WriteResults(double decayratio);
	
};


#endif
