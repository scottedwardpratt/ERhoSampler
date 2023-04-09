#include "msu_commonutils/constants.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_sampler/sampler.h"
#include "msu_erhosampler/erhosampler.h"

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	Crandy *randy=new Crandy(time(NULL));
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/parameters.txt");
	int NEVENTS_TOT=parmap.getI("NEVENTS_TOT",10);
	Eigen::MatrixXd chi4test(4,4);
	Eigen::VectorXd EQtot(7),EQTarget(7);
	for(int a=0;a<7;a++){
		EQtot(a)=0.0;
		EQTarget(a)=1.0;
	}
	
	CresList *reslist=new CresList(&parmap);
	Csampler::reslist=reslist;
	Csampler::randy=randy;
	Csampler::parmap=&parmap;
	Csampler::CALCMU=true;
	Csampler::SETMU0=false;
	Csampler::USE_POLE_MASS=true;
	Csampler::mastersampler=nullptr;
	CresInfo::randy=randy;
	
	Chyper *hyper=new Chyper();
	CcorrVsEtaScott corrvseta;
	CcorrVsY corrvsy(&parmap);
	
	int ievent;
	double T=0.150,tau=10.0,R=5.0,deleta=0.05;
	double rhoB=0.1,rhoII=0.03;
	Csampler *sampler=new Csampler(T,0.093);
	
	CpartList *partlista=new CpartList(&parmap,reslist);
	CpartList *partlistb=new CpartList(&parmap,reslist);
	
	NMSU_ERrhoSampler::FillOutHyperBjorken(hyper,T,tau,R,deleta,rhoB,rhoII);
	hyper->sampler=sampler;
	
	sampler->CalcDensitiesMu0();
	sampler->GetNHMu0();
	sampler->GetMuNH(hyper);
	
	sampler->CalcChiSlow(hyper);
	
	for(ievent=0;ievent<NEVENTS_TOT;ievent++){
		sampler->partlist=partlista;
		
		sampler->MakeParts(hyper);
		
		partlista->SetEQWeightVec(hyper);
		//DecayParts(randy,partlista);
		NMSU_ERrhoSampler::Chi4Test(partlistb,chi4test);
		
		//IncrementQtest(partlista,EQtot,EQTarget);
		partlista->TestEQWeights(EQtot,EQTarget);
		NMSU_ERrhoSampler::Chi4Test(partlista,chi4test);
		
		sampler->partlist=partlistb;
		sampler->MakeParts(hyper);
		partlistb->SetEQWeightVec(hyper);
		//DecayParts(randy,partlistb);
		
		//IncrementQtest(partlistb,EQtot,EQTarget);
		partlistb->TestEQWeights(EQtot,EQTarget);
		NMSU_ERrhoSampler::Chi4Test(partlistb,chi4test);
		
		corrvsy.Increment(partlista,partlistb,&corrvseta);
		partlista->Clear();
		partlistb->Clear();
		
		if(10*(ievent+1)%NEVENTS_TOT==0){
			printf("finished %d percent\n",((ievent+1)*100)/NEVENTS_TOT);
		}
	}
	
	chi4test=chi4test/(2*hyper->udotdOmega*NEVENTS_TOT);
	printf("------- chi4test --------\n");
	cout << chi4test << endl;
	printf("-------------------------\n");
		
	EQtot=EQtot/(2*NEVENTS_TOT);
	cout << "EQtot=\n" << EQtot << endl;
	
	corrvsy.WriteResults();
	
	delete partlista;
	delete partlistb;
	
	return 0;
}
