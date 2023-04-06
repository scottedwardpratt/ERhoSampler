#include "msu_ERhoSampler/master.h"
#include "msu_commonutils/constants.h"
#include "../software/src/ERhoSubs.cc"
#include "../software/src/CorrSubs.cc"
using namespace std;

// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	Eigen::VectorXd EQTarget(7),EQtot(7);
	int i;
	for(i=0;i<7;i++){
		EQTarget(i)=1.0;
		EQtot(i)=0.0;
	}
	double T=0.150,tau=10.0,R=5.0,eta=0.0,deleta=0.05;
	double epsilon=0.3,rhoB=0.1,rhoII=0.03;
	Csampler *sampler;
	long long int npartstot=0,ievent;
	int nparts;
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/parameters.txt");
	CmasterSampler ms(&parmap);
	CpartList *partlista=new CpartList(&parmap,ms.reslist);
	CpartList *partlistb=new CpartList(&parmap,ms.reslist);
	
	ms.randy->reset(time(NULL));
	ms.MakeDummyHyper(1);
	Chyper *hyper=*(ms.hyperlist.begin());
	FillOutHyperBjorken(hyper,T,tau,R,eta,deleta,epsilon,rhoB,rhoII);
	//hyper->Print();
	
	sampler=ms.ChooseSampler(hyper);
	sampler->GetNHMu0();
	sampler->GetMuNH(hyper);
	hyper->sampler=sampler;
	
	sampler->CalcChi(hyper);
	hyper->Print();
	cout << hyper->chi << endl;
	
	for(ievent=0;ievent<ms.NEVENTS_TOT;ievent++){
		ms.partlist=partlista;
		npartsa=ms.MakeEvent();
		ms.partlist->SetEQWeight(hyper,EQTarget);
		ms.partlist->IncrementEQTot(EQtot);
		
		
		ms.partlist=partlistb;
		npartsb=ms.MakeEvent();
		ms.partlist->SetEQWeight(hyper,EQTarget);
		ms.partlist->IncrementEQTot(EQtot);
		
		if(10*(ievent+1)%ms.NEVENTS_TOT==0){
			printf("finished %lld percent\n",((ievent+1)*100)/ms.NEVENTS_TOT);
		}
	}
	
	
	
	ms.ClearHyperList();
	
	for(i=0;i<7;i++){
		printf("%8.6f\n",EQtot(i)/double(npartstot));
	}
	
	delete partlist;
	return 0;
}
