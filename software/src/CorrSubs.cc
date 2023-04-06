double GetPairWeight(Cpart *part1,Cpart *part2,Eigen::MatrixXd &CorrMatrix);
	
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
			if(a!=b){
				if(a==1 || b==1 || a==2 || b==2)
					corr(a,b)=0.0;
				else if (a==3 || b==3){
					corr(a,b)*=eta;
					if(a==3)
						corr(a,b)*=-1;
				}
			}
		}
	}
}

CcorrVsY::CcorrVsY(){
	int a,b,iy;
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
void IncrementCorrVsY(CpartList *partlista,CpartList *partlistb,CcorrVsEta *corrvseta,CcorrVsY *corrvsy,Crandy *randy){
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
	
	for(iy=0;iy<corrvsy->NY;iy++){
		DelY=(iy+0.5)*corrvsy->DY;
	
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
		
		
				// For TESTING
				ua[0]=cosh(-ya0);
				ua[3]=sinh(-ya0);
				ub[0]=cosh(DelY-yb0);
				ub[3]=sinh(DelY-yb0);
				Misc::Boost(ua,parta->p,pa);
				Misc::Boost(ub,partb->p,pb);
			
				/*double DYTest=atanh(pa[3]/pa[0])-atanh(pb[3]/pb[0]);
				printf("DelY=%g =? %g, ya=%g, yb=%g \n",DYTest,DelY,atanh(pa[3]/pa[0]),atanh(pb[3]/pb[0]));
				Misc::Pause();*/
			
				resinfoa=parta->resinfo;
				resinfob=partb->resinfo;
				Qa[0]=pa[0];
				Qa[1]=pa[1];
				Qa[2]=pa[2];
				Qa[3]=pa[3];
				Qa[4]=resinfoa->baryon;
				Qa[5]=resinfoa->q[0]-resinfoa->q[1];
				Qa[6]=resinfoa->strange;
				Qb[0]=pb[0];
				Qb[1]=pb[1];
				Qb[2]=pb[2];
				Qb[3]=pb[3];
				Qb[4]=resinfob->baryon;
				Qb[5]=resinfob->q[0]-resinfob->q[1];
				Qb[6]=resinfob->strange;
		
				corrvsy->denom[iy]+=1.0;
	
				for(a=0;a<7;a++){
					for(b=0;b<7;b++){
						corrvsy->corr[iy](a,b)+=Qa[a]*Qb[b]*cweight;
					}
				}
			}
			//
		}
	}
}

void CcorrVsY::WriteResults(double decayratio){
	int a,b,iy,sign;
	double DelY,C,sum=0.0;
	FILE *fptr;
	char filename[120];
	for(a=0;a<7;a++){
		for(b=0;b<7;b++){
			sign=1;
			if(a>0 && a<4 && (b==0 || b>=4))
				sign=-1;
			if(b>0 && b<4 && (a==0 || a>=4))
				sign=-1;
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
