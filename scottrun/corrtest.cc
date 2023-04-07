#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"

#include "../software/src/ERhoSubs.cc"
#include "../software/src/CorrSubs.cc"


using namespace std;

// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	Chyper *hyper=new Chyper();
	CcorrVsEta corrvseta;
	CcorrVsY corrvsy;
	long long int Ndecay=0,Noriginal=0;
	Crandy *randy=new Crandy(time(NULL));
	int ievent;
	double T=0.150,tau=10.0,R=5.0,deleta=0.05;
	double rhoB=0.1,rhoII=0.03;
	Csampler *sampler;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/parameters.txt");
	CmasterSampler ms(&parmap);
	CpartList *partlista=new CpartList(&parmap,ms.reslist);
	CpartList *partlistb=new CpartList(&parmap,ms.reslist);
	
	ms.randy->reset(time(NULL));
	ms.hyperlist.push_back(hyper);
	
	FillOutHyperBjorken(hyper,T,tau,R,deleta,rhoB,rhoII);
	//ms.MakeDummyHyper(1);
	//Chyper *hyper=*(ms.hyperlist.begin());
	
	cout << hyper << " =? " << *(ms.hyperlist.begin()) << endl;
	
	//FillOutHyperBjorken(hyper,T,tau,R,deleta,rhoB,rhoII);
	
	sampler=ms.ChooseSampler(hyper);
	hyper->sampler=sampler;
	
	cout << sampler << " =? " << hyper->sampler << endl;
	
	sampler->CalcDensitiesMu0();
	sampler->GetNHMu0();
	sampler->GetMuNH(hyper);

	//sampler->CalcChi(hyper);
	sampler->GetEpsilonRhoChi(hyper->muB,hyper->muII,hyper->muS,hyper->epsilon,hyper->rhoB,hyper->rhoII,hyper->rhoS,hyper->chi);
	printf("------ chi --------\n");
	cout << hyper->chi << endl;
	
	for(ievent=0;ievent<ms.NEVENTS_TOT;ievent++){
		ms.partlist=partlista;
		ms.MakeEvent();
		Noriginal+=partlista->nparts;
		ms.partlist->SetEQWeightVec(hyper);
		DecayParts(randy,partlista);
		Ndecay+=partlista->nparts;
		
		ms.partlist=partlistb;
		ms.MakeEvent();
		Noriginal+=partlistb->nparts;
		ms.partlist->SetEQWeightVec(hyper);
		DecayParts(randy,partlistb);
		Ndecay+=partlistb->nparts;
		
		IncrementCorrVsY(partlista,partlistb,&corrvseta,&corrvsy,randy);
		
		if(10*(ievent+1)%ms.NEVENTS_TOT==0){
			printf("finished %lld percent\n",((ievent+1)*100)/ms.NEVENTS_TOT);
		}
	}
	
	corrvsy.WriteResults(double(Ndecay)/double(Noriginal));
	
	ms.ClearHyperList();

	delete partlista;
	delete partlistb;
	
	return 0;
}
