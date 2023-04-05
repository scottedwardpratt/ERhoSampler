double GetPairWeight(Cpart *part1,Cpart *part2,Eigen::MatrixXd &CorrMatrix);
	
class CcorrVsEta{
public:
	double DETA;
	void GetCorrVsEta(double eta,Eigen::MatrixXd &corr);
};

void CcorrVsEta::GetCorrVsEta(double eta,Eigen::MatrixXd &corr){
	int a,b;
	for(a=0;a<7;a++){
		for(b=0;b<7;b++){
			corr(a,b)=exp(-eta*eta/2.0);
		}
	}
}

// Dely=yb-ya
void IncrementCorrVsY(CpartList *partlista,CpartList *partlistb,double DelY,CcorrVsEta *corrvseta,Eigen::MatrixXd &CorrMatrix,Eigen::Matrix Xd &DenomMatrix,Crandy *randy){
	Eigen::MatrixXd corrmatrix(7,7);
	Eigen::VectorXd Qa(7),Qb(7);
	FourVector pa,pb,ua,ub;
	CresInfo *resinfoa,*resinfob;
	long long int imc,NMC=1000000;
	int npartsa,npartsb,ia,ib,a,b;
	double ya0,yb0,DelEta,cweight;
	Cpart *parta,*partb;
	npartsa=partlista->nparts;
	npartsb=partlistb->nparts;
	ua[1]=ua[2]=ub[1]=ub[2]=0.0;
	
	for(imc=0;imc<NMC;imc++){
		ia=lrint(floor(randy->ran()*npartsa));
		ib=lrint(floor(randy->ran()*npartsb));
		parta=&(partlista->partvec[ia]);
		partb=&(partlistb->partvec[ib]);
		ya0=parta->GetRapidity();
		yb0=partb->GetRapidity();
		DelEta=DelY+ya0-yb0;
		corrvseta->GetCorrVsEta(DelEta,corrmatrix);
		cweight=GetPairWeight(parta,partb,corrmatrix);
		
		ua[0]=cosh(-ya0);
		ua[3]=sinh(-ya0);
		ub[0]=cosh(DelY-yb0);
		ub[3]=sinh(DelY-yb0);
		Misc::Boost(ua,parta->p,pa);
		Misc::Boost(ub,partb->p,pb);
		resinfoa=parta->resinfo;
		resinfob=partb->resinfo;
		Qa[0]=parta->p[0];
		Qa[1]=parta->p[1];
		Qa[2]=parta->p[2];
		Qa[3]=parta->p[3];
		Qa[4]=resinfoa->baryon;
		Qa[5]=resinfoa->q[0]-resinfoa->q[1];
		Qa[6]=resinfoa->strange;
		Qb[0]=partb->p[0];
		Qb[1]=partb->p[1];
		Qb[2]=partb->p[2];
		Qb[3]=partb->p[3];
		Qb[4]=resinfob->baryon;
		Qb[5]=resinfob->q[0]-resinfob->q[1];
		Qb[6]=resinfob->strange;
		
		for(a=0;a<7;a++){
			for(b=0;b<7;b++){
				CorrMatrix(a,b)+=Qa[a]*Qb[b]*cweight;
				DenomMatrix(a,b)+=Qa[a]*Qb[b];
			}
		}
	}
}

void PrintCorrVsY(double DelY,CcorrVsEta *corrvseta,Eigen::MatrixXd &CorrMatrix,Eigen::Matrix Xd &DenomMatrix){
	int a,b;
	
	
}
