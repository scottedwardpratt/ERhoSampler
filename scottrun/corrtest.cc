#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"
#include "msu_erhosampler/erhosampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	int ievent;
	double T=0.150,tau=11.0,A=100.0,deleta=0.01;
	double rhoB=(8.0/11.0)*0.16,rhoQ;
	rhoQ=0.4*rhoB;
	
	Crandy *randy=new Crandy(time(NULL));
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/parameters.txt");
	int NEVENTS_TOT=parmap.getI("NEVENTS_TOT",10);
	
	CresList *reslist=new CresList(&parmap);
	CresInfo::randy=randy;
	
	Chyper *hyper=new Chyper();
	//CcorrVsEtaScott corrvseta;
	//CcorrVsEtaOleh corrvseta;
	CcorrVsEtaOlehAlt corrvseta;
	CcorrVsY corrvsy(&parmap);
	
	//Eigen::MatrixXd corrmatrix(7,7);
	//int a=4,b=4;
	// double eta;
	//for(eta=-1+0.5*0.1;eta<1.0;eta+=0.1){
		//corrvseta.GetCorrVsEta(eta,corrmatrix);
		//printf("%6.3f %8.5f\n",eta,corrmatrix(a,b));
	//}
	
	Csampler *sampler=new Csampler(T,0.093,&parmap,reslist,randy);

	CpartList *partlista=new CpartList(&parmap,reslist);
	CpartList *partlistb=new CpartList(&parmap,reslist);
	
	NMSU_ERrhoSampler::FillOutHyperBjorken(hyper,T,tau,A,deleta,rhoB,rhoQ);
	hyper->sampler=sampler;
	
	sampler->CalcDensitiesMu0();
	sampler->GetNHMu0();
	sampler->GetMuNH(hyper);
	sampler->CalcChi4BQS(hyper);
	cout << "chi4BQS\n";
	cout << hyper->chi4BQS << endl;
	
	corrvseta.TestSumRules(hyper);
	exit(1);
	
	for(ievent=0;ievent<NEVENTS_TOT;ievent++){
		sampler->partlist=partlista;
		sampler->MakeParts(hyper);
		
		
		partlista->SetEQWeightVec(hyper);
		NMSU_ERrhoSampler::DecayParts(randy,partlista);
		
		
		sampler->partlist=partlistb;
		sampler->MakeParts(hyper);
		
		partlistb->SetEQWeightVec(hyper);
		NMSU_ERrhoSampler::DecayParts(randy,partlistb);
		
		corrvsy.Increment(partlista,partlistb,&corrvseta);
		partlista->Clear();
		partlistb->Clear();
		
		if(10*(ievent+1)%NEVENTS_TOT==0){
			printf("finished %d percent\n",((ievent+1)*100)/NEVENTS_TOT);
		}
	}
	
	corrvsy.WriteResults();
	
	delete partlista;
	delete partlistb;
	
	return 0;
}
