#include<bits/stdc++.h>
#include<iostream>
#include<algorithm>
#include<cmath>
#include <iomanip>
#include<fstream>
using namespace std;
/**************************GLOBAL VARIABLE************************/
const int n=11;                                               //// pls change the array size if it is not 5      INPUT
double A[n]={0},B[n]={0},S[n]={0}; 
double Xw[n]={0}, Xe[n]={0}, X[n]={0}, dxw[n]={0}, dxe[n]={0}, dx[n]={0};
double Pp[n]={0.0}, h_inf[n]={0.0}, rho[n]={10.0}, T_inf[n]={0.0};                     //INPUT
double CAE[n]={0}, CAEO[n]={0}, CAW[n]={0}, CAWO[n]={0}, CAQ[n]={0}, CAQO[n]={0}, CAC[n]={0}, CACO[n]={0}, AT[n]={0}, ATO[n]={0}, CAP[n]={0}, CAPO[n]={0};
double Area[n]={0.0},Area_Face[n+1]={0.0},Th_Cond[n]={0.0},Th_Cond_Face[n+1]={0.0};
double Temp[n+1]={0},Temp_old[n+1]={0},Time=0;
double psi=1.0;
double mf;
double tl;
double tr;
double q_l= -10000.0;
double q_r= 0.00;
double kind_r;
double kind_l;
double Cp=750.0;
double dt=1.0;
double q[n]={0.0};
double l = 0.1;		// width of slab            INPUT
double R=2e-3;
double kind_si;
double s = l / (2.0 * (n - 1));

/******************************************************************/

template < typename T, size_t n > void
print_array (T const (&arr)[n]){
  for (size_t i = 0; i < n; i++)
    {
      std::cout<<setprecision(15)<< arr[i] << ',';
    }
}
void Gen_Face (){				//deefining Area_Face and Th_Cond_Face     INPUT--------------------INTERFACES AREA AND THERMAL CONDUCTIVITY-----------------
  for (int i = 0; i <= n; i++)
    {
      Area_Face[i] = 1.00;
      Th_Cond_Face[i] = 20.0;
    }
}
void Gen (){			//defining Area and Th_Cond     INPUT----------------------AT NODES--------------------------------------
  for (int i = 0; i < n; i++)
    {
    if (i==0){
      Area[i] =1.00;  //22*R*R/7;
      Th_Cond[i] = 20.0;
      Pp[i] = 0.00;
      //q[i]=0.25e6;

    }
    else if(i==n-1){
      Area[i] =1.00; //(22/7)*R*R;
      Th_Cond[i] = 20.0;
      Pp[i] = 0.00;
      //q[i]=0.25e6;
 
    }
    else{
        Area[i] =1.00;//(22/7)*R*R;
      Th_Cond[i] = 20.0;
      Pp[i] = 0.00;
     // h_inf[i] = 12.50;
     // T_inf[i] = 20.00;
      //q[i]=0.25e6;
    }
    }
}


//----------------------------------------------------------GRID GENERATION------------------------------------------------------------------------
void gridGeneration ( )	//Grid Generation
{
  Xw[0] = 0;			// xw[]
  for (int i = 1; i < n; i++)
    {
      if (i == 1)
	{
	  Xw[i] = Xw[i - 1] + s;
	}
      else
	{
	  Xw[i] = Xw[i - 1] + 2 * s;
	}
    }

  Xe[0] = s;			// xe[]
  for (int i = 1; i < n; i++)
    {
      if (i == n - 1)
	{
	  Xe[i] = Xe[i - 1] + s;
	}
      else{
	  Xe[i] = Xe[i - 1] + 2 * s;
	}
    }

  X[0] = 0;			// x[]
  for (int i = 1; i < n; i++){
      X[i] = X[i - 1] + 2 * s;}

  for (int i = 0; i < n; i++)	{             // dxe dxw dx
    if (i == 0){
      dxw[i] = 0;
      dx[i] = s;
      dxe[i] = 2 * s;}
    
    else if (i == n - 1)
      {
        dxw[i] = 2 * s;
        dx[i] = s;
        dxe[i] = 0;
      }
          else
      {
        dxw[i] = 2 * s;
        dxe[i] = 2 * s;
        dx[i] = 2 * s;
      }
    }
  cout << "\n";

  cout << " xw= ";
  print_array (Xw);
  cout << "\n";
  cout << "\n";
  cout << " xe= ";
  print_array (Xe);
  cout << "\n";
  cout << "\n";
  cout << " x= ";
  print_array (X);
  cout << "\n";
  cout << "\n";
  cout << "dxw= ";
  print_array (dxw);
  cout << "\n";
  cout << "\n";
  cout << "dxe= ";
  print_array (dxe);
  cout << "\n";
  cout << "\n";
  cout << "dx= ";
  print_array (dx);
  cout << "\n";
  

}

//------------------------------------------------------COEFFICIENT MATRIX CALCULATION-----------------------------------------------------------------
void CoefficientMatrix ()	
{
    for (int i = 1; i < n - 1; i++){
        CAE[i] = (psi * Area_Face[i + 1] * Th_Cond_Face[i + 1]) / dxe[i];
	      CAEO[i] = ((1 - psi) * Area_Face[i + 1] * Th_Cond_Face[i + 1]) / dxe[i];
	      CAW[i] = (psi * Area_Face[i] * Th_Cond_Face[i]) / dxw[i];
        CAWO[i] = ((1 - psi) * Area_Face[i] * Th_Cond_Face[i]) / dxw[i];
        CAQ[i] = (psi * (Area[i] * q[i] * dx[i]));
        CAQO[i] = (1 - psi) * (Area[i] * q[i] * dx[i]);
        CAC[i] = (psi * (Pp[i] * h_inf[i] * dx[i]));
        CACO[i] = (1 - psi) * (Pp[i] * h_inf[i] * dx[i]);
        AT[i] = (mf * (rho[i] * Area[i] * Cp * dx[i]) / dt);
        ATO[i] = (mf * (rho[i] * Area[i] * Cp * dx[i]) / dt);
    }
    for (int i = 1; i < n - 1; i++){
      CAP[i] = (AT[i] + CAW[i] + CAE[i] + CAC[i]);
      CAPO[i] = (ATO[i] - (CAWO[i] + CAEO[i] + CACO[i]));                           // for intermediate nodes only---------------------------------
    }
  /* cout << "\n";
  cout << " CAW= ";
  print_array (CAW);
  cout << "\n";
  cout << "\n";
  cout << " CAE= ";
  print_array (CAE);
  cout << "\n";
  cout << "\n";
  cout << " CAP= ";
  print_array (CAP);
  cout << "\n";
  cout << " CAC= ";
  print_array (CAC);
  cout << "\n";
  cout << "\n";
  cout << " CAQ= ";
  print_array(CAQ); */
 
}

void RightBoundary(){               //right BCs---coefficient-------------------------------------------------
    if (kind_r==1){
        CAE[n - 1] = 0;		//BCs       INPUT
        CAW[n - 1] = 0;
        CAP[n - 1] = 1;		//BCs       INPUT
        AT[n-1] = (mf * (rho[n-1] * Area[n-1] * Cp * dx[n-1]) / dt);
        ATO[n-1] = (mf * (rho[n-1] * Area[n-1] * Cp * dx[n-1]) / dt);        
    }
    if (kind_r==2.0){
        CAEO[n-1] = 0;		//BCs       INPUT
        CAE[n-1] = 0;
        CAW[n-1] = (psi * Area_Face[n-1] * Th_Cond_Face[n-1]) / dxw[n-1]; 
        CAC[n-1] = (psi * (Pp[n-1] * h_inf[n-1] * dx[n-1]));
        CAQ[n-1] = (psi * (Area[n-1] * q[n-1] * dx[n-1]));  
        AT[n-1] = (mf * (rho[n-1] * Area[n-1] * Cp * dx[n-1]) / dt);
        CAP[n-1] = (AT[n-1] + CAW[n-1] + CAE[n-1] + CAC[n-1]);
        ATO[n-1] = (mf * (rho[n-1] * Area[n-1] * Cp * dx[n-1]) / dt);
        CAPO[n-1] = (ATO[n-1] - (CAWO[n-1] + CAEO[n-1] + CACO[n-1]));
    }
    if (kind_r==3.0){
        CAEO[n-1] = 0;		//BCs       INPUT
        CAE[n-1] = 0;
        CAW[n-1] = (psi * Area_Face[n-1] * Th_Cond_Face[n-1]) / dxw[n-1]; 
        CAC[n-1] = (psi * (Pp[n-1] * h_inf[n-1] * dx[n-1]));
        CAQ[n-1] = (psi * (Area[n-1] * q[n-1] * dx[n-1]));  
        AT[n-1] = (mf * (rho[n-1] * Area[n-1] * Cp * dx[n-1]) / dt);
        CAP[n-1] = (AT[n-1] + CAW[n-1] + CAE[n-1] + CAC[n-1]) + Area_Face[n]*h_inf[n-1] ;
        ATO[n-1] = (mf * (rho[n-1] * Area[n-1] * Cp * dx[n-1]) / dt);
        CAPO[n-1] = (ATO[n-1] - (CAWO[n-1] + CAEO[n-1] + CACO[n-1]));
    }
}
void LeftBoundary(){               //left BCs------coefficients-------------------------------------------------
    if (kind_l==1.0){
        CAE[0] = 0;		//BCs       INPUT
        CAW[0] = 0;
        CAP[0] = 1;		//BCs       INPUT
        AT[0] = (mf * (rho[0] * Area[0] * Cp * dx[0]) / dt);
        ATO[0] = (mf * (rho[0] * Area[0] * Cp * dx[0]) / dt);
        CAC[0] = (psi * (Pp[0] * h_inf[0] * dx[0]));
        CAQ[0] = (psi * (Area[0] * q[0] * dx[0]));
        CAPO[0] = (ATO[0] - (CAWO[0] + CAEO[0] + CACO[0]));
    }
    if(kind_l==2.0){
        CAWO[0] = 0;		//BCs       INPUT
        CAW[0] = 0;
        CAE[0] = (psi * Area_Face[0] * Th_Cond_Face[0]) / dxe[0]; 
        CAC[0] = (psi * (Pp[0] * h_inf[0] * dx[0]));
        CAQ[0] = (psi * (Area[0] * q[0] * dx[0]));
        AT[0] = (mf * (rho[0] * Area[0] * Cp * dx[0]) / dt);
        CAP[0] = (AT[0] + CAW[0] + CAE[0] + CAC[0]);
        ATO[0] = (mf * (rho[0] * Area[0] * Cp * dx[0]) / dt);
        CAPO[0] = (ATO[0] - (CAWO[0] + CAEO[0] + CACO[0]));
    }
    if(kind_l==3.0){
        CAWO[0] = 0;		//BCs       INPUT
        CAW[0] = 0;
        CAE[0] = (psi * Area_Face[0] * Th_Cond_Face[0]) / dxe[0]; 
        CAC[0] = (psi * (Pp[0] * h_inf[0] * dx[0]));
        CAQ[0] = (psi * (Area[0] * q[0] * dx[0]));
        AT[0] = (mf * (rho[0] * Area[0] * Cp * dx[0]) / dt);
        CAP[0] = (AT[0] + CAW[0] + CAE[0] + CAC[0]) + Area_Face[0]*h_inf[0];
        ATO[0] = (mf * (rho[0] * Area[0] * Cp * dx[0]) / dt);
        CAPO[0] = (ATO[0] - (CAWO[0] + CAEO[0] + CACO[0]));
    }
}
//---------------------------------------------------SOURCE VECTOR CALCULATION----------------------------------------------------------------------
void sourceVectorCalculation (){				//sourceVectorCalculation(double t0,double tn)

//S[0]=tl;
//S[n-1]=S[n-1]=((CAPO[n-1] * Temp_old[n-1]) + CAQ[n-1] + CAC[n-1] * T_inf[n-1]) + psi*Area_Face[n]*h_inf[n-1]*T_inf[n-1] ;
//S[0]=((CAPO[0] * Temp_old[0]) + CAQ[0] + CAC[0] * T_inf[0]) +  (psi * Area_Face[0] * q_l );

  for (int i = 1; i < n-1 ; i++){
    S[i] =((CAPO[i] * Temp_old[i]) + CAQ[i] + CAC[i] * T_inf[i]);
    }               // sourceVector for all nodes
    // boundary conditions---------------------------------------------------------------------------------------------------------------------------
    if (q_l!=0){
        S[0]=((CAPO[0] * Temp_old[0]) + CAQ[0] + CAC[0] * T_inf[0])+  (psi * Area_Face[0] * q_l + (1 - psi) * Area_Face[0] * q_l);
    } 
    else if(tl!=0.0 && kind_si!=1){
        S[0]=tl;
    }  
    else if(h_inf[0]!=0.0){
    S[0]=((CAPO[0] * Temp_old[0]) + CAQ[0] + CAC[0] * T_inf[0])+ psi*Area_Face[0]*h_inf[0]*T_inf[0] + (1-psi)*Area_Face[0]*h_inf[0]*T_inf[0];
    } 
    else if(kind_si==1.0){
      S[0]=((CAPO[0] * Temp_old[0]) + CAQ[0] + CAC[0] * T_inf[0]);    //if left boundary is insulated.....tl is not constant && ql=0
    }              
    if(q_r!=0){
        S[n-1]=((CAPO[n-1] * Temp_old[n-1]) + CAQ[n-1] + CAC[n-1] * T_inf[n-1]) - (psi * Area_Face[n-1] * q_r - (1 - psi) * Area_Face[n-1] * q_r);
    }                                                                   
    else if(tr!=0 && kind_si!=2.0){
        S[n-1]=tr;
    }                                                      
                                                                                           
    else if(h_inf[n-1]!=0.0){
        S[n-1]=((CAPO[n-1] * Temp_old[n-1]) + CAQ[n-1] + CAC[n-1] * T_inf[n-1])+ psi*Area_Face[n-1]*h_inf[n-1]*T_inf[n-1] + (1-psi)*Area_Face[n-1]*h_inf[n-1]*T_inf[n-1];
    } 
    else if (kind_si==2.0){
      S[n-1]=((CAPO[n-1] * Temp_old[n-1]) + CAQ[n-1] + (CAC[n-1] * T_inf[n-1]));          // if right boundary insulated -----tr is not mantained and qr=0
    }                                                                                     
 
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------insulation BC-------------------------------------------------------------------------------------
   
    
}
//---------------------------------MARTIX INVERSION()----------------------------------------------------------------------------------
void matrixInversion (){				
// forward substitution
  for (int i = 0; i < n; i++){
      if (i == 0)	{
        A[i] = CAE[i] / CAP[i];
	      B[i] = S[i] / CAP[i];
	    }
      else{
        A[i] = CAE[i] / (CAP[i] - CAW[i] * A[i - 1]);
        B[i] = (S[i] + CAW[i] * B[i - 1]) / (CAP[i] - CAW[i] * A[i - 1]);
      }
  }
  cout<<"S= ";
  print_array(S);
  cout<<"\n";
  cout<<"ATO= ";
  print_array(ATO);
  cout<<"\n";
  cout<<"A= ";
  print_array(A);
  cout<<"\n";
  cout<<"B= ";
  print_array(B);
  cout<<"\n";
  cout<<"CAP= ";
  print_array(CAP);
  cout<<"\n";
  cout<<"\n";
  cout<<"CAPO= ";
  print_array(CAPO);
  cout<<"\n";
  
  cout<<"AT= ";
  print_array(AT);
  cout<<"\n";
  cout<<"CAC= ";
  print_array(CAC);
  cout<<"\n";
  cout<<"CAW= ";
  print_array(CAW);
  cout<<"\n";
  cout<<"CAQ= ";
  print_array(CAQ);
  cout<<"\n";

  
  for (int i = n - 1; i >= 0; i--)
    {
      Temp[i] = A[i] * Temp[i + 1] + B[i];
    }
  
  cout << "the T[]= " ;
  for (int i=0;i<n;i++){
    cout<<setprecision(8)<<Temp[i]<<",";
  
    }
    cout<<endl;
    //print_array (Temp);
  /*for(int i=0;i<n;i++)
  cout<<Temp[i]<<",";*/
  cout<<endl;

}
//-----------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------TransientSolver --------------------------------------------------------------------------------
void TransientSolver()
{
  mf=1.0;
  tl=150.000;
  tr=75.000;
  for(int i=1;i<n-1;i++)
  {
    Temp_old[i]=tl+(tr-tl)*X[i]/l;
    //Temp[i]=tl+(tr-tl)*X[i]/l;
  }
 
  Temp_old[0]=tl;
  Temp_old[n-1]=tr;
  //Temp[0]=tl;
  //Temp[n-1]=tr;
  cout<<"\n";
  cout<<"\n";
  cout<<"before time loop starts ";
  cout<<endl;
  print_array(Temp_old);
  cout<<"\n";
  h_inf[n-1]=100.00; 
  T_inf[n-1]=175.00;



  for(int i=1;i<=3000;i++)                                                           // specify time level
  {
    tl=0.0;
    tr=0.0;
    Time+=dt;
    cout<<"n="<<i<<" Time="<<Time;
    cout<<"\n";

    kind_l=2.0;       // type of left BC
    kind_r=3.0;       // type of right BC
    kind_si=0.0;      // type -where is insulation-left=1 or right=2

    RightBoundary();
    LeftBoundary();
    CoefficientMatrix();
    
    sourceVectorCalculation();
    matrixInversion();
    
    for(int j=0;j<n;j++){
    Temp_old[j]=Temp[j];}

  }
  
  
}
//-----------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------SteadySolver --------------------------------------------------------------------------------
void SteadySolver()
{
  mf=0;
   kind_l=1.0;       // type of left BC
    kind_r=1.0;       // type of right BC
    kind_si=0.0;      // type -where is insulation-left=1 or right=2
  RightBoundary();
  LeftBoundary();
  CoefficientMatrix ();
  sourceVectorCalculation ();	
  matrixInversion ();
  for(int j=0;j<n;j++)
    Temp_old[j]=Temp[j];
}
//------------------------------------MAIN PROGRAM------------------------------------------------------------------------------
int main ()
{
  
  h_inf[0]=0.0;
 
  //T_inf[0]=20.0;
 


  for(int i=0;i<n;i++)
    rho[i]=2000.0;
  
  gridGeneration ();
  Gen_Face ();
  Gen ();
  
  
  //SteadySolver();
  TransientSolver();

  
  
  //cout<<"Temp_old= "<<Temp_old[1];
 
  //cout<<"mf= "<<mf;
  ofstream graval;
   graval.open("GraphValue.csv",ios::out);
   graval<<"Xw,Xe,X,dxw,dxe,dx,CAE,CAEO, CAW, CAWO, CAQ, CAQO, CAC, CACO, AT, ATO, CAP, CAPO, S, T"<<endl;

   for(int i=0;i<n;i++)
   {
        graval<<Xw[i]<<","<<Xe[i]<<","<<X[i]<<","<<dxw[i]<<","<<dxe[i]<<","<<dx[i]<<","<<CAE[i]<<","<<CAEO[i]<<","<< CAW[i]<<","<< CAWO[i]<<","<< CAQ[i]<<","<< CAQO[i]<<","<< CAC[i]<<","<< CACO[i]<<","<< AT[i]<<","<< ATO[i]<<","<< CAP[i]<<","<< CAPO[i]<<","<< S[i]<<","<< Temp[i]<<endl;
   }
return 0;
}