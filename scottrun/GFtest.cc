#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>

using namespace std;


// This makes a dummy hyper-element then creates particles and tests yield and energy of created partilces of specific pid

int main(){
	int ievent;
	double eta,etamin=-5.0,etamax=5.0,deta=0.02,dtau=0.02,tau0,tau00=1.0,tauf=11.0;
	int ieta,itau0,NETA,NTAU;
	double sumtestBB,sumtestBE,sumtestEE,sumtestBP,sumtestEP,sumtestPP,sumtestPE,csheta,ssheta;
	vector<vector<double>> GF_BB,GF_BE,GF_EE,GF_BP,GF_PP,GF_EP,GF_PE;
	NETA=lrint((etamax-etamin)/deta)+1;
	NTAU=lrint((tauf-tau00)/dtau);
	FILE *fptr;
	
	GF_BB.resize(NTAU);
	GF_BE.resize(NTAU);
	GF_EE.resize(NTAU);
	GF_BP.resize(NTAU);
	GF_PP.resize(NTAU);
	GF_EP.resize(NTAU);
	GF_PE.resize(NTAU);
	for(itau0=0;itau0<NTAU;itau0++){
		GF_BB[itau0].resize(NETA);
		GF_BE[itau0].resize(NETA);
		GF_EE[itau0].resize(NETA);
		GF_BP[itau0].resize(NETA);
		GF_PP[itau0].resize(NETA);
		GF_EP[itau0].resize(NETA);
		GF_PE[itau0].resize(NETA);
		for(ieta=0;ieta<NETA;ieta++){
			GF_BB[itau0][ieta]=0.0;
			GF_BE[itau0][ieta]=0.0;
			GF_EE[itau0][ieta]=0.0;
			GF_BP[itau0][ieta]=0.0;
			GF_PP[itau0][ieta]=0.0;
			GF_EP[itau0][ieta]=0.0;
			GF_PE[itau0][ieta]=0.0;
		}
	}

	// B-B
	fptr=fopen("OlehCFs/Green-B-B.dat","r");
	for(itau0=0;itau0<NTAU;itau0++){
		tau0=tau00+dtau*itau0;
		for(ieta=0;ieta<NETA;ieta++){
			eta=etamin+ieta*deta;
			fscanf(fptr,"%lf",&GF_BB[itau0][ieta]);
			if(feof(fptr)){
				printf("file ended: itau0=%d, ieta=%d\n",itau0,ieta);
				exit(1);
			}
			GF_BB[itau0][ieta]*=tauf;
			eta=etamin+ieta*deta;
		}
	}
	fclose(fptr);
	
	printf("check a\n");
	// B-E
	fptr=fopen("OlehCFs/Green-B-E.dat","r");
	for(itau0=0;itau0<NTAU;itau0++){
		tau0=tau00+dtau*itau0;
		for(ieta=0;ieta<NETA;ieta++){
			eta=etamin+ieta*deta;
			fscanf(fptr,"%lf",&GF_BE[itau0][ieta]);
			if(feof(fptr)){
				printf("file ended: itau0=%d, ieta=%d\n",itau0,ieta);
				exit(1);
			}
			GF_BE[itau0][ieta]*=tauf;
			eta=etamin+ieta*deta;
		}
	}
	fclose(fptr);
	
	printf("check b\n");
	// E-E
	fptr=fopen("OlehCFs/Green-E-E.dat","r");
	for(itau0=0;itau0<NTAU;itau0++){
		tau0=tau00+dtau*itau0;
		for(ieta=0;ieta<NETA;ieta++){
			eta=etamin+ieta*deta;
			fscanf(fptr,"%lf",&GF_EE[itau0][ieta]);
			if(feof(fptr)){
				printf("file ended: itau0=%d, ieta=%d\n",itau0,ieta);
				exit(1);
			}
			GF_EE[itau0][ieta]*=tauf;
			eta=etamin+ieta*deta;
		}
	}
	fclose(fptr);
	
	printf("check c\n");
	// B-P
	fptr=fopen("OlehCFs/Green-B-Uz.dat","r");
	for(itau0=0;itau0<NTAU;itau0++){
		tau0=tau00+dtau*itau0;
		for(ieta=0;ieta<NETA;ieta++){
			eta=etamin+ieta*deta;
			fscanf(fptr,"%lf",&GF_BP[itau0][ieta]);
			if(feof(fptr)){
				printf("file ended: itau0=%d, ieta=%d\n",itau0,ieta);
				exit(1);
			}
			GF_BP[itau0][ieta]*=tauf;
			eta=etamin+ieta*deta;
		}
	}
	fclose(fptr);
	
	printf("check d\n");
	// P-P
	fptr=fopen("OlehCFs/Green-U-Uz.dat","r");
	for(itau0=0;itau0<NTAU;itau0++){
		tau0=tau00+dtau*itau0;
		for(ieta=0;ieta<NETA;ieta++){
			eta=etamin+ieta*deta;
			fscanf(fptr,"%lf",&GF_PP[itau0][ieta]);
			if(feof(fptr)){
				printf("file ended: itau0=%d, ieta=%d\n",itau0,ieta);
				exit(1);
			}
			GF_PP[itau0][ieta]*=tauf;
			eta=etamin+ieta*deta;
		}
	}
	fclose(fptr);
	
	printf("check e\n");
	// E-P
	fptr=fopen("OlehCFs/Green-U-E.dat","r");
	for(itau0=0;itau0<NTAU;itau0++){
		tau0=tau00+dtau*itau0;
		for(ieta=0;ieta<NETA;ieta++){
			eta=etamin+ieta*deta;
			fscanf(fptr,"%lf",&GF_PE[itau0][ieta]);
			if(feof(fptr)){
				printf("file ended: itau0=%d, ieta=%d\n",itau0,ieta);
				exit(1);
			}
			GF_PE[itau0][ieta]*=tauf;
			eta=etamin+ieta*deta;
		}
	}
	fclose(fptr);
	
	printf("check f\n");
	// E-P
	fptr=fopen("OlehCFs/Green-E-Uz.dat","r");
	for(itau0=0;itau0<NTAU;itau0++){
		tau0=tau00+dtau*itau0;
		for(ieta=0;ieta<NETA;ieta++){
			eta=etamin+ieta*deta;
			fscanf(fptr,"%lf",&GF_EP[itau0][ieta]);
			if(feof(fptr)){
				printf("file ended: itau0=%d, ieta=%d\n",itau0,ieta);
				exit(1);
			}
			GF_EP[itau0][ieta]*=tauf;
			eta=etamin+ieta*deta;
		}
	}
	fclose(fptr);
	
	
	// Do sum tests
	for(itau0=0;itau0<NTAU;itau0++){
		tau0=tau00+dtau*itau0;
		sumtestBB=sumtestBE=sumtestEE=sumtestBP=sumtestPP=sumtestEP=sumtestEP=0.0;
		for(ieta=0;ieta<NETA;ieta++){
			eta=etamin+ieta*deta;
			csheta=cosh(eta);
			ssheta=sinh(eta);
			sumtestBB+=deta*GF_BB[itau0][ieta];
			sumtestEE+=deta*(GF_EE[itau0][ieta]*csheta+GF_EP[itau0][ieta]*ssheta);
			sumtestEP+=deta*(GF_EP[itau0][ieta]*csheta+GF_PP[itau0][ieta]*ssheta);
			sumtestBE+=deta*(GF_BE[itau0][ieta]*csheta+GF_BP[itau0][ieta]*ssheta);
			sumtestBP+=deta*(GF_BP[itau0][ieta]*csheta+GF_BE[itau0][ieta]*ssheta);
			sumtestPP+=deta*(GF_PP[itau0][ieta]*csheta+GF_PE[itau0][ieta]*ssheta);
		}
		printf("tau0=%5.3f  sumtestBB=%6.3f, sumtestEE=%6.3f, sumtestEP=%6.3f, sumtestBE=%6.3f, sumtestBP=%6.3f, sumtestPP=%6.3f\n",tau0,sumtestBB,sumtestEE,sumtestEP,sumtestBE,sumtestBP,sumtestPP);
	}
	
	return 0;
}
