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
	int imother,idaughter;
	int nbodies,ibody;
	double mtot;
	vector<Cpart> daughterpart(10);
	vector<double> masses(11);
	vector<FourVector> pdaughter(10);
	Cpart *motherpart,*dpart;
	array<CresInfo *,5> daughterresinfo;
	imother=0;
	while(imother<partlist->nparts){
		motherpart=&(partlist->partvec[imother]);
		if(motherpart->resinfo->decay){
			motherpart->resinfo->DecayGetResInfoPtr(nbodies,daughterresinfo);
			Decay.nbodies=nbodies;
			mtot=0.0;
			masses[0]=motherpart->resinfo->mass;
			for(ibody=0;ibody<nbodies;ibody++){
				masses[ibody+1]=daughterresinfo[ibody]->mass;
				mtot+=masses[ibody+1];
			}
			// Increase mother mass if not large enough
			if(mtot>masses[0]){
				masses[0]=mtot+0.01;
			}
			motherpart->msquared=masses[0]*masses[0];
			motherpart->Setp0();
			mtot=0.0;
			for(ibody=0;ibody<nbodies;ibody++)
				mtot+=masses[ibody+1];
			if(mtot>masses[0]){
				printf("This cannot be, m[0]=%g, m[1]=%g, m[2]=%g, mtot=%g, nbodies=%d\n",masses[0],masses[1],masses[2],mtot,nbodies);
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
					dpart->p=pdaughter[ibody];
					dpart->SetMsquared();
				}
				else{
					partlist->AddPart(daughterresinfo[ibody]->pid,pdaughter[ibody],motherpart->r);
					idaughter=partlist->nparts-1;
					dpart=&(partlist->partvec[idaughter]);
					dpart->EQWeightVec=motherpart->EQWeightVec;
				}
			}
		}
		if(!motherpart->resinfo->decay)
			imother++;
	}
}
