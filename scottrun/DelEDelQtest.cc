#include "msu_ERhoSampler/master.h"
#include "msu_commonutils/constants.h"
using namespace std;

// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

void FillOutHyperBjorken(Chyper *hyper,double T,double tau,double R,double eta,double deleta,
double epsilon,double rhoB,double rhoII){
	double V0=PI*R*R*tau*deleta;
	hyper->T0=T;
	hyper->sigma=0.093;
	hyper->rhoB=rhoB;
	hyper->rhoS=0.0;
	hyper->rhoII=rhoII;
	hyper->epsilon=epsilon;
	hyper->muB=hyper->muS=hyper->muII=0.0;
	hyper->u[1]=hyper->u[2]=hyper->u[3]=sinh(eta);
	hyper->u[0]=sqrt(1.0+hyper->u[1]*hyper->u[1]+hyper->u[2]*hyper->u[2]+hyper->u[3]*hyper->u[3]);
	hyper->r[1]=hyper->r[2]=hyper->r[3]=0.0;
	hyper->r[0]=10.0;
	for(int alpha=0;alpha<4;alpha++)
		hyper->dOmega[alpha]=V0*hyper->u[alpha];//*2.0*ms.randy->ran();
	hyper->udotdOmega=hyper->u[0]*hyper->dOmega[0]-hyper->u[1]*hyper->dOmega[1]
		-hyper->u[2]*hyper->dOmega[2]-hyper->u[3]*hyper->dOmega[3];
	for(int alpha=0;alpha<4;alpha++)
		for(int beta=0;beta<4;beta++)
			hyper->pitilde[alpha][beta]=0.0;
}

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
	parmap.ReadParsFromFile("parameters.txt");
	CmasterSampler ms(&parmap);
	CpartList *partlist=new CpartList(&parmap,ms.reslist);
	ms.partlist=partlist;
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
		// Perform events
		nparts=ms.MakeEvent();
		npartstot+=nparts;
		ms.partlist->SetEQWeight(hyper,EQTarget);
		ms.partlist->IncrementEQTot(EQtot);
		
		if(10*(ievent+1)%ms.NEVENTS_TOT==0){
			printf("finished %lld percent\n",((ievent+1)*100)/ms.NEVENTS_TOT);
			double nparts_target=hyper->nhadrons*hyper->udotdOmega*double(ievent);
			printf("npartstot=%lld =? %g, ratio=%g\n",npartstot,nparts_target,npartstot/nparts_target);
		}
	}
	ms.ClearHyperList();
	
	for(i=0;i<7;i++){
		printf("%8.6f\n",EQtot(i)/double(npartstot));
	}
	
	delete partlist;
	return 0;
}
