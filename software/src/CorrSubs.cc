#include "msu_erhosampler/erhosampler.h"

class CcorrVsEta{
public:
	double DETA;
	void GetCorrVsEta(double eta,Eigen::MatrixXd &corr);
};

class CcorrVsY{
public:
	const double DY=0.1;
	const int NY=100;
	vector<Eigen::MatrixXd> corr;
	vector<double> denom;
	void WriteResults(double decayratio);
	CcorrVsY();
};

void CcorrVsEta::GetCorrVsEta(double eta,Eigen::MatrixXd &corr){
	int a,b;
	double Z=1.0/sqrt(2.0*PI);
	for(a=0;a<7;a++){
		for(b=0;b<7;b++){
			corr(a,b)=Z*exp(-eta*eta/2.0);
		}
	}
	//corr(0,3)=eta*Z*exp(-eta*eta/2.0);
	//corr(3,0)=-eta*Z*exp(-eta*eta/2.0);
}

CcorrVsY::CcorrVsY(CparameterMap *parmap_set){
	parmap=parmap_set;
	int a,b,iy;
	DY=parmap->getD("MSU_ERHOSAMPLER_DY",0.1);
	NY=parmap->getI("MSU_ERHOSAMPLER_DY",100);
	corr.resize(NY);
	denom.resize(NY);
	for(iy=0;iy<NY;iy++){
		corr[iy].resize(7,7);
		for(a=0;a<7;a++)
			for(b=0;b<7;b++)
				corr[iy](a,b)=0.0;
	}
}

// Dely=yb-ya
void CcorrVsY::Increment(CpartList *partlista,CpartList *partlistb,CcorrVsEta *corrvseta){
	Eigen::MatrixXd corrmatrix(7,7);
	Eigen::VectorXd Qa(7),Qb(7);
	FourVector pa,pb,ua,ub;
	CresInfo *resinfoa,*resinfob;
	//long long int imc,NMC=lrint(partlista->nparts*partlistb->nparts);
	int npartsa,npartsb,a,b,iy,ia,ib;
	double ya0,yb0,DelEta,cweight,DelY;
	Cpart *parta,*partb;
	npartsa=partlista->nparts;
	npartsb=partlistb->nparts;
	ua[1]=ua[2]=ub[1]=ub[2]=0.0;
	
	for(iy=0;iy<NY;iy++){
		DelY=(iy+0.5)*DY;
	
		//for(imc=0;imc<NMC;imc++){
		//ia=lrint(floor(randy->ran()*npartsa));
		//ib=lrint(floor(randy->ran()*npartsb));
		for(ia=0;ia<npartsa;ia++){
			for(ib=0;ib<npartsb;ib++){
				parta=&(partlista->partvec[ia]);
				partb=&(partlistb->partvec[ib]);
				ya0=parta->GetRapidity();
				yb0=partb->GetRapidity();
				DelEta=DelY+ya0-yb0;
				corrvseta->GetCorrVsEta(DelEta,corrmatrix);
				cweight=GetPairWeight(parta,partb,corrmatrix);
		
		
				/* For TESTING
				ua[0]=cosh(-ya0);
				ua[3]=sinh(-ya0);
				ub[0]=cosh(DelY-yb0);
				ub[3]=sinh(DelY-yb0);
				Misc::Boost(ua,parta->p,pa);
				Misc::Boost(ub,partb->p,pb);*/
			
				resinfoa=parta->resinfo;
				resinfob=partb->resinfo;
				//Qa[0]=pa[0];
				//Qa[1]=pa[1];
				//Qa[2]=pa[2];
				//Qa[3]=pa[3];
				Qa[0]=parta->p[0];
								Qa[1]=parta->p[1];
												Qa[2]=parta->p[2];
																Qa[3]=parta->p[3];
				Qa[4]=resinfoa->baryon;
				Qa[5]=resinfoa->q[0]-resinfoa->q[1];
				Qa[6]=resinfoa->strange;
				//Qb[0]=pb[0];
				//Qb[1]=pb[1];
				//Qb[2]=pb[2];
				//Qb[3]=pb[3];
				Qb[0]=partb->p[0];
								Qb[1]=partb->p[1];
												Qb[2]=partb->p[2];
																Qb[3]=partb->p[3];
				Qb[4]=resinfob->baryon;
				Qb[5]=resinfob->q[0]-resinfob->q[1];
				Qb[6]=resinfob->strange;
		
				denom[iy]+=1.0;
	
				for(a=0;a<7;a++){
					for(b=0;b<7;b++){
						corr[iy](a,b)+=Qa[a]*Qb[b]*cweight;
					}
				}
			}
			//
		}
	}
}


double CcorrVsY::GetPairWeight(Cpart *part1,Cpart *part2,Eigen::MatrixXd &CorrMatrix){
	int i,j;
	double weight=0.0;
	for(i=0;i<7;i++){
		for (j=0;j<7;j++){
			weight+=part1->EQWeightVec[i]*CorrMatrix(i,j)*part2->EQWeightVec[j];
		}
	}
	return weight;
}

void CcorrVsY::WriteResults(double decayratio){
	int a,b,iy,sign;
	double DelY,C,sum=0.0;
	FILE *fptr;
	char filename[120];
	for(a=0;a<7;a++){
		for(b=0;b<7;b++){
			sign=1;
			/*if(a>0 && a<4 && (b==0 || b>=4))
				sign=-1;
			if(b>0 && b<4 && (a==0 || a>=4))
				sign=-1;*/
			sum=0.0;
			snprintf(filename,120,"corr_results/corr_%d_%d.txt",a,b);
			fptr=fopen(filename,"w");
			for(iy=0;iy<NY;iy++){
				DelY=DY*(iy+0.5);
				C=decayratio*decayratio*(corr[iy](a,b)+sign*corr[iy](b,a))/denom[iy];
				fprintf(fptr,"%6.3f %g\n",DelY,C);
				sum+=C*DY;
			}
			printf("sum[%d][%d]=%g\n",a,b,sum);
			fclose(fptr);
		}
	}
}

void IncrementQtest(CpartList *partlist,Eigen::VectorXd &Qtot,Eigen::VectorXd &EQTarget){
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
		Q[5]=resinfo->q[0]-resinfo->q[1];
		Q[6]=resinfo->strange;
		weight=0.0;
		for(a=0;a<7;a++){
			weight+=part->EQWeightVec[a]*EQTarget[a];
		}
		for(a=0;a<7;a++)
			Qtot[a]+=Q[a]*weight;
	}
}

void Chi4Test(CpartList *partlist,Eigen::MatrixXd &chitest){
	int a,b,ia,nparts=partlist->nparts;
	Eigen::VectorXd Q(4);
	Cpart *part;
	CresInfo *resinfo;
	for(ia=0;ia<nparts;ia++){
		part=&(partlist->partvec[ia]);
		resinfo=part->resinfo;
		Q[0]=part->p[0];
		Q[1]=resinfo->baryon;
		Q[2]=resinfo->q[0]-resinfo->q[1];
		Q[3]=resinfo->strange;
		for(a=0;a<4;a++){
			for(b=0;b<4;b++){
				chitest(a,b)+=Q(a)*Q(b);
			}
		}
	}
}
