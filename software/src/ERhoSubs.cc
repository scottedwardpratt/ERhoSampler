#include "msu_erhosampler/erhosampler.h"

void NMSU_ERrhoSampler::FillOutHyperBjorken(Chyper *hyper,double T,double tau,double A,double deleta,double rhoB,double rhoQ){
	double V0=A*tau*deleta;
	hyper->T0=T;
	hyper->sigma=0.093;
	hyper->rhoB=rhoB;
	hyper->rhoS=0.0;
	hyper->rhoQ=rhoQ;
	hyper->rhoII=2.0*rhoQ-rhoB;
	hyper->epsilon=0.3;
	hyper->muB=hyper->muS=hyper->muII=0.0;
	hyper->u[1]=hyper->u[2]=hyper->u[3]=0.0;
	hyper->u[0]=sqrt(1.0+hyper->u[1]*hyper->u[1]+hyper->u[2]*hyper->u[2]+hyper->u[3]*hyper->u[3]);
	hyper->r[1]=hyper->r[2]=hyper->r[3]=0.0;
	hyper->r[0]=tau;
	hyper->tau=tau;
	hyper->sampler=nullptr;
	for(int alpha=0;alpha<4;alpha++)
		hyper->dOmega[alpha]=V0*hyper->u[alpha];//*2.0*ms.randy->ran();
	hyper->udotdOmega=hyper->u[0]*hyper->dOmega[0]-hyper->u[1]*hyper->dOmega[1]
		-hyper->u[2]*hyper->dOmega[2]-hyper->u[3]*hyper->dOmega[3];
	for(int alpha=0;alpha<4;alpha++)
		for(int beta=0;beta<4;beta++)
			hyper->pitilde[alpha][beta]=0.0;
}

void NMSU_ERrhoSampler::DecayParts(Crandy *randyset,CpartList *partlist){
	CDecay_NBody Decay(randyset);
	int imother,idaughter,alpha,a;
	int nbodies,ibody,ntry;
	Eigen::VectorXd EQWeightVecTmp(7);
	double mtot;
	char message[CLog::CHARLENGTH];
	FourVector pmother,pafter,pprime;
	vector<Cpart> daughterpart;
	vector<double> masses;
	vector<FourVector> pdaughter;
	pdaughter.resize(10);
	daughterpart.resize(10);
	masses.resize(11);
	int qbefore[3],qafter[3];
	FourVector u;
	Cpart *motherpart,*dpart;
	array<CresInfo *,5> daughterresinfo;
	imother=0;
	while(imother<partlist->nparts){
		motherpart=&(partlist->partvec[imother]);
		EQWeightVecTmp=motherpart->EQWeightVec;
		pmother=motherpart->p;
		for(alpha=0;alpha<3;alpha++){
			pafter[alpha]=0.0;
			qafter[alpha]=0;
			qbefore[alpha]=motherpart->resinfo->q[alpha];
		}
		if(motherpart->resinfo->decay){
			motherpart->Setp0();
			masses[0]=motherpart->resinfo->mass;
			ntry=0;
			do{
				motherpart->resinfo->DecayGetResInfoPtr(nbodies,daughterresinfo);
				
				mtot=0.0;
				for(ibody=1;ibody<=nbodies;ibody++){
					masses[ibody]=daughterresinfo[ibody-1]->mass;
					mtot+=masses[ibody];
				}
				ntry+=1;
				if(ntry>100){
					motherpart->resinfo->Print();
					CLog::Info("mothermass="+to_string(masses[0])+", minmass="+to_string(motherpart->resinfo->minmass)+"\n");
					CLog::Fatal("cannot find channel\n");
				}
			}while(mtot>masses[0]);
			
			for(alpha=0;alpha<4;alpha++)
				u[alpha]=motherpart->p[alpha]/masses[0];
			//
			Decay.SetMasses(nbodies,masses);
			Decay.GenerateMomenta(pdaughter);
			for(ibody=0;ibody<nbodies;ibody++){
				if(ibody==0){
					dpart=motherpart;
					dpart->pid=daughterresinfo[ibody]->pid;
					dpart->resinfo=daughterresinfo[ibody];
				}
				else{
					partlist->AddPart(dpart->pid=daughterresinfo[ibody]->pid,pdaughter[ibody],motherpart->r);
					idaughter=partlist->nparts-1;
					dpart=&(partlist->partvec[idaughter]);
					dpart->EQWeightVec=motherpart->EQWeightVec;
				}
				dpart->EQWeightVec=EQWeightVecTmp;
				dpart->p=pdaughter[ibody];
				dpart->SetMsquared();
				Misc::Boost(u,dpart->p,pprime);
				for(alpha=0;alpha<4;alpha++)
					dpart->p[alpha]=pprime[alpha];
				for(a=0;a<3;a++){
					qafter[a]+=dpart->resinfo->q[a];
				}
				for(alpha=0;alpha<4;alpha++){
					pafter[alpha]+=dpart->p[alpha];
				}
			}
			for(a=0;a<3;a++){
				if(qbefore[a]!=qafter[a]){
					snprintf(message,CLog::CHARLENGTH,"qbefore=(%d,%d,%d),qafter=(%d,%d,%d)\n",qbefore[0],qbefore[1],qbefore[2],
					qafter[0],qafter[1],qafter[2]);
					CLog::Fatal(message);
				}
			}
		}
		if(!motherpart->resinfo->decay)
			imother++;
	}
}

void NMSU_ERrhoSampler::IncrementEQtest(CpartList *partlist,Eigen::VectorXd &EQtot,Eigen::VectorXd &EQTarget){
	int a,ia,nparts=partlist->nparts;
	Eigen::VectorXd Q(7);
	Cpart *part;
	CresInfo *resinfo;
	double weight;
	for(ia=0;ia<nparts;ia++){
		part=&(partlist->partvec[ia]);
		resinfo=part->resinfo;
		Q[0]=part->p[0];
		Q[1]=part->p[1];
		Q[2]=part->p[2];
		Q[3]=part->p[3];
		Q[4]=resinfo->baryon;
		Q[5]=resinfo->charge;
		Q[6]=resinfo->strange;
		weight=0.0;
		for(a=0;a<7;a++){
			weight+=part->EQWeightVec[a]*EQTarget[a];
		}
		for(a=0;a<7;a++)
			EQtot[a]+=Q[a]*weight;
	}
}

void NMSU_ERrhoSampler::Chi4Test(CpartList *partlist,Eigen::MatrixXd &chitest){
	int a,b,ia,nparts=partlist->nparts;
	Eigen::VectorXd Q(4);
	Cpart *part;
	CresInfo *resinfo;
	for(ia=0;ia<nparts;ia++){
		part=&(partlist->partvec[ia]);
		resinfo=part->resinfo;
		Q[0]=part->p[0];
		Q[1]=resinfo->baryon;
		Q[2]=resinfo->charge;
		Q[3]=resinfo->strange;
		for(a=0;a<4;a++){
			for(b=0;b<4;b++){
				chitest(a,b)+=Q(a)*Q(b);
			}
		}
	}
}
