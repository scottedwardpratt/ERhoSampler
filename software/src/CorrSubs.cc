#include "msu_erhosampler/erhosampler.h"

using namespace NMSU_ERrhoSampler;

CcorrVsY::CcorrVsY(CparameterMap *parmap_set){
	parmap=parmap_set;
	int a,b,iy;
	denom=0.0;
	DY=parmap->getD("MSU_ERHOSAMPLER_DY",0.1);
	NY=parmap->getI("MSU_ERHOSAMPLER_NY",200);
	corr.resize(NY);
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
	int npartsa,npartsb,a,b,iy,ia,ib;
	double ya0,yb0,DelEta,cweight,DelY;
	Cpart *parta,*partb;
	npartsa=partlista->nparts;
	npartsb=partlistb->nparts;
	denom+=1;
	ua[1]=ua[2]=ub[1]=ub[2]=0.0;
	
	for(iy=0;iy<NY;iy++){
		DelY=(-(NY/2)+iy+0.5)*DY;
		for(ia=0;ia<npartsa;ia++){
			for(ib=0;ib<npartsb;ib++){
				parta=&(partlista->partvec[ia]);
				partb=&(partlistb->partvec[ib]);
				ya0=parta->GetRapidity();
				yb0=partb->GetRapidity();
				DelEta=DelY+ya0-yb0;
				corrvseta->GetCorrVsEta(DelEta,corrmatrix);
				for(int a=0;a<7;a++)
					for(int b=0;b<7;b++)
						if(a!=0 || b!=0)
							corrmatrix(a,b)=0.0;
				cweight=GetPairWeight(parta,partb,corrmatrix);
			
				resinfoa=parta->resinfo;
				resinfob=partb->resinfo;
				
				pa[3]=0.0;
				pa[0]=sqrt(parta->p[0]*parta->p[0]-parta->p[3]*parta->p[3]);
				Qa[0]=pa[0];
				Qa[1]=pa[1];
				Qa[2]=pa[2];
				Qa[3]=pa[3];
				//Qa[0]=parta->p[0];
				//Qa[1]=parta->p[1];
				//Qa[2]=parta->p[2];
				//Qa[3]=parta->p[3];
				pb[3]=0.0;
				pb[0]=sqrt(partb->p[0]*partb->p[0]-partb->p[3]*partb->p[3]);
				Qa[4]=resinfoa->baryon;
				Qa[5]=resinfoa->charge;
				Qa[6]=resinfoa->strange;
				Qb[0]=pb[0];
				Qb[1]=pb[1];
				Qb[2]=pb[2];
				Qb[3]=pb[3];
				//Qb[0]=partb->p[0];
				//Qb[1]=partb->p[1];
				//Qb[2]=partb->p[2];
				//Qb[3]=partb->p[3];
				Qb[4]=resinfob->baryon;
				Qb[5]=resinfob->charge;
				Qb[6]=resinfob->strange;
	
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

void CcorrVsY::WriteResults(){
	int a,b,iy;
	double DelY,C,sum=0.0;
	FILE *fptr;
	char filename[120];
	for(a=0;a<7;a++){
		for(b=0;b<7;b++){
			sum=0.0;
			snprintf(filename,120,"hydrocorr_results/corr_%d_%d.txt",a,b);
			fptr=fopen(filename,"w");
			for(iy=0;iy<NY;iy++){
				DelY=(-(NY/2)+iy+0.5)*DY;
				C=corr[iy](a,b)/denom;
				fprintf(fptr,"%6.3f %g\n",DelY,C);
				sum+=C*DY;
			}
			printf("sum[%d][%d]=%g\n",a,b,sum);
			fclose(fptr);
		}
	}
}
