// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

void FillOutHyperBjorken(Chyper *hyper,double T,double tau,double R,double deleta,double rhoB,double rhoII){
	double V0=PI*R*R*tau*deleta;
	hyper->T0=T;
	hyper->sigma=0.093;
	hyper->rhoB=rhoB;
	hyper->rhoS=0.0;
	hyper->rhoII=rhoII;
	hyper->epsilon=0.3;
	hyper->muB=hyper->muS=hyper->muII=0.0;
	hyper->u[1]=hyper->u[2]=hyper->u[3]=0.0;
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

void DecayParts(Crandy *randyset,CpartList *partlist){
	CDecay_NBody Decay(randyset);
	int imother,idaughter,alpha;
	int nbodies,ibody;
	vector<Cpart> daughterpart;
	vector<double> masses;
	vector<FourVector> pdaughter;
	pdaughter.resize(10);
	daughterpart.resize(10);
	masses.resize(11);
	FourVector u;
	Cpart *motherpart,*dpart;
	array<CresInfo *,5> daughterresinfo;
	imother=0;
	while(imother<partlist->nparts){
		motherpart=&(partlist->partvec[imother]);
		if(motherpart->resinfo->decay){
			motherpart->resinfo->DecayGetResInfoPtr(nbodies,daughterresinfo);
			motherpart->Setp0();
			masses[0]=motherpart->resinfo->mass;
			for(ibody=0;ibody<nbodies;ibody++){
				masses[ibody+1]=daughterresinfo[ibody]->mass;
			}
			
			double mtot=0.0;
			for(ibody=1;ibody<=nbodies;ibody++)
				mtot+=masses[ibody+1];
			if(mtot>masses[0]){
				printf("YIKES!!!!\n");
				exit(1);
			}
			
			for(alpha=0;alpha<4;alpha++)
				u[alpha]=motherpart->p[alpha]/masses[0];
			double u2=u[0]*u[0]-u[1]*u[1]-u[2]*u[2]-u[3]*u[3];
			if(fabs(u2-1.0)>0.0001){
				printf("1=?%g, mtot=%g, mothermass=%g\\\n",u2,mtot,masses[0]);
				printf("p_mother=(%g,%g,%g,%g)\n",motherpart->p[0],motherpart->p[1],motherpart->p[2],motherpart->p[3]);
				exit(1);
			}
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
				dpart->p=pdaughter[ibody];
				dpart->SetMsquared();
				Misc::Boost(u,dpart->p,dpart->p);
			}
		}
		if(!motherpart->resinfo->decay)
			imother++;
	}
}

double GetPairWeight(Cpart *part1,Cpart *part2,Eigen::MatrixXd &Corr){
	int i,j;
	double weight=0.0;
	for(i=0;i<7;i++){
		for (j=0;j<7;j++){
			weight+=part1->EQWeightVec[i]*Corr(i,j)*part2->EQWeightVec[j];
		}
	}
	return weight;
}
