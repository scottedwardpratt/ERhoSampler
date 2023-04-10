#include "msu_erhosampler/erhosampler.h"

using namespace NMSU_ERrhoSampler;

void CcorrVsEtaScott::GetCorrVsEta(double eta,Eigen::MatrixXd &corr){
	int a,b;
	double Z=1.0/sqrt(2.0*PI);
	for(a=0;a<7;a++){
		for(b=0;b<7;b++){
			corr(a,b)=0.0;
			if((a==0 ||a>3) && (b==0||b>3)){
				corr(a,b)=Z*exp(-eta*eta/2.0);
			}
			else if(a==b){
				corr(a,b)=Z*exp(-eta*eta/2.0);
			}
			else if(a==3 && (b==0||b>3)){
				corr(a,b)=Z*eta*exp(-eta*eta/2.0);
			}
			else if(b==3 && (a==0||a>3)){
				corr(a,b)=-Z*eta*exp(-eta*eta/2.0);
			}
		}
	}
	//corr(0,3)=eta*Z*exp(-eta*eta/2.0);
	//corr(3,0)=-eta*Z*exp(-eta*eta/2.0);
}

void CcorrVsEtaOleh::GetCorrVsEta(double eta,Eigen::MatrixXd &corr){
	int a,b;
	double Z=1.0/sqrt(2.0*PI);
	for(a=0;a<7;a++){
		for(b=0;b<7;b++){
			corr(a,b)=0.0;
			if((a==0 ||a>3) && (b==0||b>3)){
				corr(a,b)=Z*exp(-eta*eta/2.0);
			}
			else if(a==b){
				corr(a,b)=Z*exp(-eta*eta/2.0);
			}
			else if(a==3 && (b==0||b>3)){
				corr(a,b)=Z*eta*exp(-eta*eta/2.0);
			}
			else if(b==3 && (a==0||a>3)){
				corr(a,b)=-Z*eta*exp(-eta*eta/2.0);
			}
		}
	}
	//corr(0,3)=eta*Z*exp(-eta*eta/2.0);
	//corr(3,0)=-eta*Z*exp(-eta*eta/2.0);
}