
#include<iostream>
#include<algorithm>
#include<cmath>
#include<fstream>
using namespace std;
/**************************GLOBAL VARIABLE************************/
const int n=5;
//double A[n]={0},B[n]={0},T[n+1]={0},S[n]={0}; 
double xw[n]={0}, xe[n]={0}, x[n]={0}, dxw[n]={0}, dxe[n]={0}, dx[n]={0};
const int m=5;      			// no of nodes         
double ys[m]={0}, yN[m]={0}, y[m]={0}, dys[m]={0}, dyN[m]={0}, dy[m]={0};
//double Pp[n]={0.0}, h_inf[n]={0.0}, T_inf[n]={0.0};                     //INPUT
//double CAE[n]={0}, CAEO[n]={0}, CAW[n]={0}, CAWO[n]={0}, CAQ[n]={0}, CAQO[n]={0}, CAC[n]={0}, CACO[n]={0}, AT[n]={0}, ATO[n]={0}, CAP[n]={0}, CAPO[n]={0};
//double Area[n]={0.0},Area_Face[n+1]={0.0},Th_Cond[n]={0.0},Th_Cond_Face[n+1]={0.0};

double lx= 0.1;		// width of slab            INPUT
double sx = lx / (2.0 * (n - 1));
double ly = 0.1;		// width of slab            INPUT
double sy = ly / (2.0 * (m - 1));

double psi=1.0;                            //implicit 
double mf=0.0;                             // steady state 

double tl= 0.0;                           // temperature at left 
double tm= 20.0;                           // temperature at right
double tn= 0.0;                           // temperature at top 
double ts= 20.0;                           // temperature at bottom

double q_l= -10000.0;
double q_r= 0.0;

double kind_r=1.0 ;                         // type of boundary at right
double kind_l=2.0;                          // type of boundary at left

//double q[n];//=400000.00;                               // heat generation
//double rho[n]={0.0};
double Cp=1000.0;
double dt=1.0;

/******************************************************************/

template < typename T, size_t n > void
print_array (T const (&arr)[n]){
  for (size_t i = 0; i < n; i++)
    {
      std::cout << arr[i] << ',';
    }
}
/*void Gen (){			//defining Area and Th_Cond     INPUT
  for (int i = 0; i < n; i++)
    {
    if (i==0){
      Area[i] = 1.00;
      Th_Cond[i] = 10.0;
      Pp[i] = 1.00;
      rho[i]=2000.0;
      q[i]=400000.00;

    }
    else if(i==n-1){
        Area[i] = 1.00;
      Th_Cond[i] = 10.0;
      Pp[i] = 1.000;
      rho[i]=2000.0;
      q[i]=400000.00;
 
    }
    else{
        Area[i] = 1.00;
      Th_Cond[i] = 10.0;
      Pp[i] = 1.00;
      h_inf[i] = 0.0;
      T_inf[i] = 0.0;
      rho[i]=2000.0;
      q[i]=400000.00;
    }
    }
}

void Gen_Face (){				//deefining Area_Face and Th_Cond_Face     INPUT
  for (int i = 0; i <= n; i++)
    {
      Area_Face[i] =1.00;
      Th_Cond_Face[i] = 10.0;
    }
}*/
//----------------------------------------------------------GRID GENERATION------------------------------------------------------------------------
void gridGeneration ()	//Grid Generation
{
  xw[0] = 0;			// xw[]
  for (int i = 1; i < n; i++)
    {
      if (i == 1)
	{
	  xw[i] = xw[i - 1] + sx;
	}
      else
	{
	  xw[i] = xw[i - 1] + 2 * sx;
	}
    }

  xe[0] = sx;			// xe[]
  for (int i = 1; i < n; i++)
    {
      if (i == n - 1)
	{
	  xe[i] = xe[i - 1] + sx;
	}
      else
	{
	  xe[i] = xe[i - 1] + 2 * sx;
	}
    }

  x[0] = 0;			// x[]
  for (int i = 1; i < n; i++)
    {
      x[i] = x[i - 1] + 2 * sx;

    }

  for (int i = 0; i < n; i++)	// dxe dxw dx
    {
      if (i == 0)
	{
	  dxw[i] = 0;
	  dx[i] = sx;
	  dxe[i] = 2 * sx;
	}
      else if (i == n - 1)
	{
	  dxw[i] = 2 * sx;
	  dx[i] = sx;
	  dxe[i] = 0;
	}
      else
	{
	  dxw[i] = 2 * sx;
	  dxe[i] = 2 * sx;
	  dx[i] = 2 * sx;
	}
    }
  cout << "\n";

  cout << " xw= ";
  print_array (xw);
  cout << "\n";
  cout << "\n";
  cout << " xe= ";
  print_array (xe);
  cout << "\n";
  cout << "\n";
  cout << " x= ";
  print_array (x);
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


  ys[0] = 0;			// xw[]
  for (int i = 1; i < m; i++)
    {
      if (i == 1)
	{
	  ys[i] = ys[i - 1] + sy;
	}
      else
	{
	  ys[i] = ys[i - 1] + 2 * sy;
	}
    }

  yN[0] = sy;			// xe[]
  for (int i = 1; i < m; i++)
    {
      if (i == m - 1)
	{
	  yN[i] = yN[i - 1] + sy;
	}
      else
	{
	  yN[i] = yN[i - 1] + 2 * sy;
	}
    }

  y[0] = 0;			// x[]
  for (int i = 1; i < m; i++)
    {
      y[i] = y[i - 1] + 2 * sy;

    }

  for (int i = 0; i < m; i++)	// dxe dxw dx
    {
      if (i == 0)
	{
	  dys[i] = 0;
	  dy[i] = sy;
	  dyN[i] = 2 * sy;
	}
      else if (i == m - 1)
	{
	  dys[i] = 2 * sy;
	  dy[i] = sy;
	  dyN[i] = 0;
	}
      else
	{
	  dys[i] = 2 * sy;
	  dyN[i] = 2 * sy;
	  dy[i] = 2 * sy;
	}
    }
  cout << "\n";

  cout << " ys= ";
  print_array (ys);
  cout << "\n";
  cout << "\n";
  cout << " yN= ";
  print_array (yN);
  cout << "\n";
  cout << "\n";
  cout << " y= ";
  print_array (y);
  cout << "\n";
  cout << "\n";
  cout << "dys= ";
  print_array (dys);
  cout << "\n";
  cout << "\n";
  cout << "dyN= ";
  print_array (dyN);
  cout << "\n";
  cout << "\n";
  cout << "dy= ";
  print_array (dy);
  cout << "\n";

}

//------------------------------------------------------COEFFICIENT MATRIx CALCULATION-----------------------------------------------------------------
/*void CoefficientMatrix ()	
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
  cout << "\n";
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

  for (int i = 0; i <= n - 1; i++){
    S[i] =((CAPO[i] * T[i]) + (CAWO[i] * T[i - 1]) + (CAEO[i] * T[i + 1]) + CAQO[i] + CAQ[i] + CAC[i] * T_inf[i] + CACO[i] * T_inf[i]);
    }               // sourceVector for all nodes
    
    // boundary conditions---------------------------------------------------------------------------------------------------------------------------
    if (q_l!=0){
        S[0]=S[0]+  (psi * Area_Face[0] * q_l + (1 - psi) * Area_Face[0] * q_l);
    }               //heat flux at left 
    if(q_r!=0){
        S[n-1]=S[n-1] - (psi * Area_Face[n-1] * q_r - (1 - psi) * Area_Face[n-1] * q_r);
    }               //heat flux at right
    if(q_l==0 && h_inf[0]==0 && tl!=0){
        S[0]=tl;
    }                                                       // temperature specified at left
    if(q_r==0 && h_inf[n-1]==0 && tm!=0){
        S[n-1]=tm;
    }                                                      // temperature specified at left
    if(h_inf[0]!=0.0){
        S[0]=S[0]+ psi*Area_Face[0]*h_inf[0]*T_inf[0] + (1-psi)*Area_Face[0]*h_inf[0]*T_inf[0];
    }                                                                                           //heat convection at left
    if(h_inf[n-1]!=0.0){
        S[n-1]=S[n-1]+ psi*Area_Face[n-1]*h_inf[n-1]*T_inf[n-1] + (1-psi)*Area_Face[n-1]*h_inf[n-1]*T_inf[n-1];
    }                                                                                          // heat convection at right
    //-------------------------------------------------------------------------------------------------------------------------------------------------
}
//---------------------------------MARTIx INVERSION()----------------------------------------------------------------------------------
void matrixInversion (){				
// forward substitution
  for (int i = 0; i < n; i++){
      if (i == 0)
	{
	  A[i] = CAE[i] / CAP[i];
	  B[i] = S[i] / CAP[i];
	}
      else
	{
	  A[i] = CAE[i] / (CAP[i] - CAW[i] * A[i - 1]);
	  B[i] = (S[i] + CAW[i] * B[i - 1]) / (CAP[i] - CAW[i] * A[i - 1]);
	}
}
  cout << "\n";
  cout << "[A]= ";
  print_array (A);
  cout << "\n";
  cout << "\n";
  cout << "[B]= ";
  print_array (B);
  cout << "\n";
//  back substitution
  for (int i = n - 1; i >= 0; i--)
    {

      T[i] = A[i] * T[i + 1] + B[i];
    }
  cout << "\n";
  cout << "the desired solution are T[n+1]=;" << endl;
  print_array (T);

}
//-----------------------------------------------------------------------------------------------------------------------------------
void heatFlux(double f){
  matrixInversion();
  double caw1,cae1,cap1,caw2,cae2,cap2,s1,s2,hFlux1,hFlux2;
  if (f==1){
  caw1=0.0;
  cae1= (psi * Area_Face[0 + 1] * Th_Cond_Face[0 + 1]) / dxe[0];
  CAC[0] = (psi * (Pp[0] * h_inf[0] * dx[0]));
  CAQ[0] = (psi * (Area[0] * q[0] * dx[0]));
  cap1=caw1 + cae1 + AT[0] + CAC[0];
  s1=cap1*T[0]-cae1*T[1];
  hFlux1=(CAQ[0]+CAC[0]*T_inf[0]-s1);
 /* cout<<'\n';
  cout<<"cap1= "<<cap1;
  cout<<'\n';
  cout<<"cae1= "<<cae1;
  cout<<'\n';
  cout<<"CAQ1= "<<CAQ[0];
  cout<<'\n';
  cout<<"s1= "<<s1;*/
 /* cout<<'\n';
  cout<<"Heat Flux at left bounadry = "<<hFlux1<<endl;
 }
  if (f==2){
  caw2=(psi * Area_Face[n-1] * Th_Cond_Face[n-1]) / dxw[n-1];
  cae2= 0.0;
  CAQ[n-1] = (psi * (Area[n-1] * q[n-1] * dx[n-1]));
  CAC[n-1] = (psi * (Pp[n-1] * h_inf[n-1] * dx[n-1]));
  cap2=caw2 + cae2 + AT[n-1] +CAC[n-1];
  s2=cap1*T[n-1]-caw1*T[n-2];
  hFlux2=(CAQ[n-1]+CAC[n-1]*T_inf[n-1]-s2);
 /* cout<<'\n';
  cout<<"cap2= "<<cap2;
  cout<<'\n';
  cout<<"caw2= "<<caw2;
  cout<<'\n';
  cout<<"CAQ2= "<<CAQ[n-1];
  cout<<'\n';
  cout<<"s2= "<<s2;*/
 /* cout<<'\n';
  cout<<"Heat Flux at right bounadry = "<<hFlux2<<endl;
  }

}*/
//------------------------------------MAIN PROGRAM------------------------------------------------------------------------------
int main ()
{
  /*h_inf[0]=0.0;
  h_inf[n-1]=0.0;
  T_inf[0]=0.0;
  T_inf[n-1]=0.0;*/

  gridGeneration ();
  //Gen_Face ();
  //Gen ();
  //RightBoundary();
  //LeftBoundary();
  //CoefficientMatrix ();
  //sourceVectorCalculation ();	
  //matrixInversion ();
  //heatFlux(2.0);            // 1 for left boundary and 2 for right boundary

 

//  print_array(S);
   ofstream graval;
 /*  graval.open("GraphValue.csv",ios::out);
   graval<<"xw,xe,x,dxw,dxe,dx,CAE,CAEO, CAW, CAWO, CAQ, CAQO, CAC, CACO, AT, ATO, CAP, CAPO, S, T"<<endl;

   for(int i=0;i<n;i++)
   {
        graval<<xw[i]<<","<<xe[i]<<","<<x[i]<<","<<dxw[i]<<","<<dxe[i]<<","<<dx[i]<<","<<CAE[i]<<","<<CAEO[i]<<","<< CAW[i]<<","<< CAWO[i]<<","<< CAQ[i]<<","<< CAQO[i]<<","<< CAC[i]<<","<< CACO[i]<<","<< AT[i]<<","<< ATO[i]<<","<< CAP[i]<<","<< CAPO[i]<<","<< S[i]<<","<< T[i]<<endl;
   }
   */  

  return 0;
}
