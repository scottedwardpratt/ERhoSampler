// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

class CcorrVsEta{
	double DETA;
	void GetCorrVsEta(double eta,Eigen::MatrixXd &corr);
};

void GetCorrVsY(CpartList *partlista,CpartList *partlistb,double DelY,CcorrVsEta *corrvseta){
	Eigen::MatrixXd &corrmatrix(7,7);
	long long int imc,NMC=1000000;
	int npartsa,npartsb,ia,ib;
	double ya0,yb0,DelEta,cweight;
	Cpart *parta,*partb;
	npartsa=partlista->nparts;
	npartsb=partlistb->nparts;
	for(imc=0;imc<NMC;imc++){
		ia=lrint(floor(randy->ran()*npartsa));
		ib=lrint(floor(randy->ran()*npartsb));
		parta=partlista->partvec[ia];
		partb=partlistb->partvec[ib];
		ya0=parta->GetRapidity();
		yb0=partb->GetRapidity();
		Deleta=DelY+yb0-ya0;
		corrvseta->GetCorrVsEta(DelETa,corrmatrix);
		cweight=GetPairWeight(parta,corrmatrix,partb);
		
	}
}
