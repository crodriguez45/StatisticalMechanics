//============================================================================
// Name        : Tarea3.cpp
// Author      : Cristian Fernando Rodríguez Cruz
// Copyright   : 
// Description : Hello World in ROOT/C++, project that works whit ROOT-libs
//============================================================================
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <unistd.h>
#include <pthread.h>
#include <time.h>
#include <fstream>
#include "TApplication.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
using namespace std;
#define NUM_THREADS 1

const int n=15;//number of rows
const int m=n;//number of columns 
const int N=n*m;
const double Kb=1;
const double J=1*Kb; 
const double Tc=2.*J/(asinh(1)*Kb);
const int NNmax=1e6;
const int NNmin=2e5;

double T=0.01;
double Beta=1./(Kb*T);
int S[n+2][m+2];
bool drawEt=false;
bool drawMt=true;
int NN=0; 
	double X[NNmax];
	double Y[NNmax];
	double ex[NNmax];
	double ey[NNmax];
	double t[NNmax];
	double Et[NNmax];
	double Mt[NNmax];

void Rand_Matrix(int a[][m+2]);
void BoundaryCond(int a[][m+2]);
void PrintMatrix(int a[][m+2]);
void Metropolis();
void DrawSp();
double DeltaE(int i, int j);
double Energy();
double Magnetization();
TCanvas *C1;
TGraphErrors *g1;
TGraph *g2;

void exec(double *En, double *mag){
	NN=0;
	srand (time(NULL));
	Rand_Matrix(S);
	BoundaryCond(S);
	//c1=new TCanvas("c1", "Lattice", 0,0,1270,720);
	*En=0;
	*mag=0;
	
	while(NN<NNmax){
		if(drawEt){
			t[NN]=NN;
			Et[NN]=Energy();
			Mt[NN]=Magnetization();
		}
		
		Metropolis();
		if(NN>NNmin){
			*En+=Energy();
			*mag+=Magnetization();
		}
		NN++;
	}
	if(drawEt){
		
		drawEt=false;
		DrawSp();
	}
	*En=*En/double(NNmax-NNmin);
	*mag=*mag/double(NNmax-NNmin);
	//delete c1;
}

int main(int argc, char **argv){

	const int nT=250;
	double Tmin=1;
	double Tmax=2.5*Tc;
	double Energia=0;
	double Magnetizacion=0;
	double TT[nT];
	double ET[nT];
	double MT[nT];
	double eTT[nT];
	double eET[nT];
	double eMT[nT];
	
	for(int i=1; i<=nT; i++){
		const int nMuestras =150;
		double X[nMuestras];
		double Y[nMuestras];
		double DesvX=0;
		double DesvY=0;
		
		T=Tmin+i*(Tmax-Tmin)/nT;
		Beta=1./(Kb*T);
		
		for(int j=0;j<nMuestras; j++){
			exec(&Energia, &Magnetizacion);
			X[j]=Energia;
			Y[j]=Magnetizacion;
		}
		
		Energia=0;
		Magnetizacion=0;
		
		for(int j=0;j<nMuestras; j++){
				Energia+=X[j]/double(nMuestras);
				Magnetizacion+=Y[j]/double(nMuestras);
		}
		
		for(int j=0;j<nMuestras; j++){
			DesvX+=pow(X[j]-Energia,2)/nMuestras;
			DesvY+=pow(Y[j]-Magnetizacion,2)/nMuestras;
		}
		
		DesvX=sqrt(DesvX);
		DesvY=sqrt(DesvY);
		
		TT[i-1]=T;
		ET[i-1]=Energia;
		MT[i-1]=Magnetizacion;
		eTT[i-1]=(Tmax-Tmin)/double(2*nT);
		eET[i-1]=DesvX;
		eMT[i-1]=DesvY;
		if( int((100*i)/nT)!=int((100*(i-1))/nT)){
			cout <<(100*i)/nT<<"% \n";
		}
	}
	
	C1= new TCanvas("C1", "Canvas", 200, 10, 700, 500);
	C1->SetFillColor(0);
	//graficar
	g1 = new TGraphErrors(nT,TT,ET,eTT,eET);
	g1->SetTitle("Energia vs Tempertura");
	g1->GetYaxis()->SetTitle("E");
	g1->GetXaxis()->SetTitle("T/J");
	g1->SetMarkerColor(4);
	g1->SetMarkerStyle(kDot);
	g1->Draw("AP");  
	C1->Update();
	C1->SaveAs("ET.pdf");
	gPad->Update();
	delete C1;
	delete g1;
	
	C1= new TCanvas("C1", "Canvas", 200, 10, 700, 500);
	C1->SetFillColor(0);
	//graficar
	g1 = new TGraphErrors(nT,TT,MT,eTT,eMT);
	g1->SetTitle("Magnetizacion vs Tempertura");
	g1->GetYaxis()->SetTitle("M");
	g1->GetXaxis()->SetTitle("T/J");
	g1->SetMarkerColor(4);
	g1->SetMarkerStyle(kDot);
	g1->Draw("AP");  
	C1->Update();
	C1->SaveAs("MT.pdf");
	gPad->Update();
	delete C1;
	delete g1;
	
return 0;
}

void DrawSp(){
	C1= new TCanvas("C1", "Canvas", 200, 10, 700, 500);
	C1->SetFillColor(0);
	g2 = new TGraph(NNmax,t,Et);
	g2->SetTitle("Energia vs iteraciones");
	g2->GetYaxis()->SetTitle("E");
	g2->GetXaxis()->SetTitle("n");
	g2->SetMarkerColor(4);
	g2->SetMarkerStyle(kDot);
	g2->Draw("AL");  
	C1->Update();
	C1->SaveAs("En.pdf");
	gPad->Update();
	delete C1;
	delete g2;
		
	C1= new TCanvas("C1", "Canvas", 200, 10, 700, 500);
	C1->SetFillColor(0);
	g2 = new TGraph(NNmax,t,Mt);
	g2->SetTitle("Magnetizacion vs n iteraciones");
	g2->GetYaxis()->SetTitle("M");
	g2->GetXaxis()->SetTitle("n");
	g2->SetMarkerColor(4);
	g2->SetMarkerStyle(kDot);
	g2->Draw("AL");  
	C1->Update();
	C1->SaveAs("Mn.pdf");
	gPad->Update();
	delete C1;
	delete g2;

	C1= new TCanvas("C1", "Canvas", 200, 10, 700, 500);
	C1->SetFillColor(0);
	TGraph2D *dt = new TGraph2D(n*m);
	dt->SetNpy(280);
	dt->SetNpx(496);
	int k=0;
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			dt->SetPoint(k,j,i,S[i][j]);
			k++;
		}
	}
	dt->SetPoint(0,0,0,-1);
	gStyle->SetPalette(kPigeon);
	dt->Draw("colz");
	C1->Update();
	C1->SaveAs("Lattice.pdf");
	delete dt;	
	delete C1;
}

double Energy(){
	double sum=0;
	for (int i=1; i<=n; i++) {
		for (int j=1; j<=m; j++) {
			sum+=S[i][j]*(S[i-1][j]+S[i][j+1]);
		}
	}
	return -J*sum/N; 
}

double Magnetization(){
	double sum=0;
	for (int i=1; i<=n; i++) {
		for (int j=1; j<=m; j++) {
			sum+=S[i][j];
		}
	}
	return abs(sum)/double(N); 
}

void Metropolis(){
	int I=rand()%n + 1;
	int J=rand()%m + 1;
	if(DeltaE(I,J)<=0){S[I][J]*=-1;}
	else{
		double p=double(rand())/double(RAND_MAX); 
		if(p<exp(-Beta*DeltaE(I,J))){
			S[I][J]*=-1;
		}
	}
	//cout<<(DeltaE(I,J))<<endl;
	if(I==1||I==n||J==1||J==m){
		BoundaryCond(S);
	}
}

double DeltaE(int i, int j){
	int DeltaS=-2*S[i][j];
	int DeltaSS=DeltaS*(S[i-1][j]+S[i+1][j]+S[i][j-1]+S[i][j+1]);
	return -J*double(DeltaSS);
}

void Rand_Matrix(int a[][m+2]){
	for(int i=1; i<=n; i++){
		for(int j=1; j<=m; j++){
			a[i][j]=2*(rand()%2)-1;
		}
	} 
}
void BoundaryCond(int a[][m+2]){
	for(int j=1; j<=m;j++){
		a[0][j]		=a[n][j];
		a[n+1][j]	=a[1][j];
	}
	for(int i=0; i<=n+1;i++){
		a[i][0]		=a[i][m];
		a[i][m+1]	=a[i][m+1];
	}
}