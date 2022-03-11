
#include<iostream>
#include<algorithm>
#include<cmath>
#include<fstream>
using namespace std;
/**************************GLOBAL VARIABLE************************/
const int n=5;      			// no of nodes                                              INPUT
double A[n]={0},B[n]={0},T[n+1]={0},S[n]={0}; 
double xw[n], xe[n], x[n], dxw[n]={0}, dxe[n]={0}, dx[n]={0};
double Pp[n]={0.0}, h_inf[n]={0.0}, T_inf[n]={0.0};                     //INPUT
double CAE[n]={0}, CAEO[n]={0}, CAW[n]={0}, CAWO[n]={0}, CAQ[n]={0}, CAQO[n]={0}, CAC[n]={0}, CACO[n]={0}, AT[n]={0}, ATO[n]={0}, CAP[n]={0}, CAPO[n]={0};
double Area[n]={0.0},Area_Face[n+1]={0.0},Th_Cond[n]={0.0},Th_Cond_Face[n+1]={0.0};



double psi=1.0;                            //implicit 
double mf=0.0;                             // steady state 

double tl= 50.0;                           // temperature at left 
double tm= 20.0;                           // temperature at right

double q_l= 0.0;
double q_r= 0.0;

double kind_r=1.0 ;                         // type of boundary at right
double kind_l=1.0;                          // type of boundary at left
double kind_si=0.0;

double q[n];//=400000.00;                               // heat generation
double rho[n]={0.0};
double Cp=1000.0;
double dt=1.0;

//----------------------------------------------------------DISCONTINUOUS GRID GENERATION------------------------------------------------------------------------
//-------------------------------------------------------------------INPUTS------------------------------------------
int const d=1;
double li[d-1];
double l[d+1]={9,6};
int M[d+1]={3,2};
double s[d+1]={0};
double L;
double A1[d+1]={10.00,10.00};
double K1[d+1]={20.00,30.00};
double AK[n+1]={0};
//--------------------------------------------------------------------------------------------------------------------

template < typename T, size_t n > void
print_array (T const (&arr)[n]){
  for (size_t i = 0; i < n; i++)
    {
      std::cout << arr[i] << ',';
    }
}
void Gen (){			//defining Area and Th_Cond     INPUT----------------------AT NODES--------------------------------------
int position =0;
  for (int i = 0; i < d+1; i++)
    for (int j=0;j<M[i];j++)
    {
        if (i==0){
        Area[position] =A1[i];  //22*R*R/7;
        Area_Face[position]=A1[i];
        Th_Cond[position] = K1[i];
        Th_Cond_Face[position]=K1[i];
        AK[position]=Area_Face[position]*Th_Cond_Face[position];
        //Pp[position] = 4*sqrt(A1[i]);
        h_inf[position] = 100.00;
        T_inf[position] = 175.00;
        //q[position]=0.25;
        position+=1;
        }
        int mid=1;
        if(i>0 && i<d){
        for(int j=0;j<M[i];j++){
            Area[position] =A1[i];  //22*R*R/7;
            Area_Face[position]=A1[i];
            Th_Cond[position] = K1[i];
            Th_Cond_Face[position]=K1[i];
            AK[position]=Area_Face[position]*Th_Cond_Face[position];
            //Pp[position] = 4*sqrt(A1[i]);
            h_inf[position] = 100.00;
            T_inf[position] = 175.00;
            
            }
            position+=1;
        }
        if(i==d){
            Area[position] =A1[i];  //22*R*R/7;
            Area_Face[position]=A1[i];
            Th_Cond[position] = K1[i];
            Th_Cond_Face[position]=K1[i];
            AK[position]=Area_Face[position]*Th_Cond_Face[position];
            //Pp[position] = 4*sqrt(A1[i]);
            h_inf[position] = 100.00;
            T_inf[position] = 175.00;
            //q[position]=0.25e6;
            position+=1;
        }
       
    }
}
void Gen_Face (){				//deefining Area_Face and Th_Cond_Face     INPUT--------------------INTERFACES AREA AND THERMAL CONDUCTIVITY-----------------
  int Flag=0;
  for(int j=M[Flag];j<=n && Flag<d+1;j+=M[Flag])
    {
      //Th_Cond_Face[j] = 2*(1/((1/Th_Cond[j-1])+(1/Th_Cond[j])));
      AK[j] = 2*(1/((1/AK[j-1])+(1/AK[j])));
      Flag++;
    }
    Th_Cond_Face[n]=Th_Cond[n-1];
    Area_Face[n]=Area[n-1];
    AK[n]=Area_Face[n]*Th_Cond_Face[n];
}

//-----------------------------------------------------GRID GENERATION------------------------------------------------
void gridGeneration()
{   //--------------------------------------------------------------------------- s CALCULATION------------------------
    //first region
    s[0]=(l[0]/(2*M[0]-1));
    //intermediate region
    int f=1;
    for(int j=1;j<d;j++){

        s[j]=((l[f]-s[j-1])/(2*M[j]-1));
        f=f+1;
    }
    //last region
    s[d]=((l[d]-s[d-1])/(2*M[d]-2));
    int index=1;
    //----------------------------------------- // xw[]----------------------------------------------------------------
    xw[0]=0.0;                       
    for (int i=0;i<d+1;i++)
    {for(int j=0;j<M[i];j++)
        {
            if (i==0 && j==0){
                xw[index]=s[i];
                index+=1;
            }
            else if(i!=0 && j==0){
                xw[index]=xw[index-1] + s[i-1] +s[i];
                index+=1;
            }
            else if (i==d && j==(M[i]-1)){
                break;
            }
            else{
                xw[index]=xw[index-1]+ 2*s[i];
                index+=1;
            }
        }
    }//------------------------------------------------------// xe[]-------------------------------------------------------
    int index2=0;              
    for (int i=0;i<d+1;i++)
    {for(int j=0;j<M[i];j++)
        {
            if (i==0 && j==0){
                xe[index2]=s[i];
                index2+=1;
            }
            else if(i!=0 && j==0){
                xe[index2]=xe[index2-1] + s[i-1] +s[i];
                index2+=1;
            }
            else if (i==d && j==(M[i]-1)){
                break;
            }
            else{
                xe[index2]=xe[index2-1]+ 2*s[i];
                index2+=1;
            }
        }
    }
    //--------------------------------------------------------- //x[]-----------------------------------------------------------
    x[0]=0.0;                  
    int index3=1;
    for (int i=0;i<d+1;i++)
    {for(int j=0;j<M[i];j++)
        {
            if (i==d && j==(M[i]-1)){
                break;
            }

            if( j==M[i]){
                x[index3]=x[index3-1] + s[i+1] + s[i];
                index3+=1;
            }
            else{
                x[index3]=x[index3-1] + 2*s[i];
                index3+=1;
            }
        }
    }//-----------------------------------------------------------//dx,dxw.dxe----------------------------------------------------
    for(int a=0;a<n;a++){       
       if (a==0){
           dx[a]=s[0];
           dxw[a]=0.0;
           dxe[a]=2*s[0];
       }
       else if(a==n-1){
           dx[a]=s[d];
           dxw[a]=2*s[d];
           dxe[a]=0.0;
       }
       else{
           dx[a]=xe[a]-xw[a];
           dxw[a]=x[a]-x[a-1];
           dxe[a]=x[a+1]-x[a];
       }
   }
  cout<<"\n";
  cout<<"xw = ";
  print_array(xw);
  
  cout<<"\n";
  cout<<"xe = ";
  print_array(xe);
  cout<<"\n";
    cout<<"x = ";
  print_array(x);
  cout<<"\n";
    cout<<"dxw = ";
  print_array(dxw);
  cout<<"\n";
    cout<<"dxe = ";
  print_array(dxe);
  cout<<"\n";
    cout<<"dx = ";
  print_array(xw);
  cout<<"\n";
   /* //-----------------------------------------------------------RESULTS------------------------------------------------------------
    cout<<setw(10)<<"xw"<<setw(10)<<"xe"<<setw(10)<<"x"<<endl;
    for (int i=0;i<n;i++)
    {cout<<"i="<<i+1<<setw(8)<<xw[i]<<setw(10)<<xe[i]<<setw(10)<<x[i]<<endl;
    }
    cout<<setw(10)<<"dxw"<<setw(10)<<"dxe"<<setw(10)<<"dx"<<endl;
    for (int i=0;i<n;i++)
    {cout<<"i="<<i+1<<setw(8)<<dxw[i]<<setw(10)<<dxe[i]<<setw(10)<<dx[i]<<endl;
    }//------------------------------------------------------------------------------------------------------------------------------- */
}
//------------------------------------------------------COEFFICIENT MATRIx CALCULATION-----------------------------------------------------------------
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

  for (int i = 1; i < n-1 ; i++){
    S[i] =((CAPO[i] * T[i]) + CAQ[i] + CAC[i] * T_inf[i]);
    }               // sourceVector for all nodes
    // boundary conditions---------------------------------------------------------------------------------------------------------------------------
    if (q_l!=0){
        S[0]=((CAPO[0] * T[0]) + CAQ[0] + CAC[0] * T_inf[0])+  (psi * Area_Face[0] * q_l + (1 - psi) * Area_Face[0] * q_l);
    } 
    else if(tl!=0.0 && kind_si!=1){
        S[0]=tl;
    }  
    else if(h_inf[0]!=0.0){
    S[0]=((CAPO[0] * T[0]) + CAQ[0] + CAC[0] * T_inf[0])+ psi*Area_Face[0]*h_inf[0]*T_inf[0] + (1-psi)*Area_Face[0]*h_inf[0]*T_inf[0];
    } 
    else if(kind_si==1.0){
      S[0]=((CAPO[0] * T[0]) + CAQ[0] + CAC[0] * T_inf[0]);    //if left boundary is insulated.....tl is not constant && ql=0
    }              
    if(q_r!=0){
        S[n-1]=((CAPO[n-1] * T[n-1]) + CAQ[n-1] + CAC[n-1] * T_inf[n-1]) - (psi * Area_Face[n-1] * q_r - (1 - psi) * Area_Face[n-1] * q_r);
    }                                                                   
    else if(tm!=0 && kind_si!=2.0){
        S[n-1]=tm;
    }                                                      
                                                                                           
    else if(h_inf[n-1]!=0.0){
        S[n-1]=((CAPO[n-1] * T[n-1]) + CAQ[n-1] + CAC[n-1] * T_inf[n-1])+ psi*Area_Face[n-1]*h_inf[n-1]*T_inf[n-1] + (1-psi)*Area_Face[n-1]*h_inf[n-1]*T_inf[n-1];
    } 
    else if (kind_si==2.0){
      S[n-1]=((CAPO[n-1] * T[n-1]) + CAQ[n-1] + (CAC[n-1] * T_inf[n-1]));          // if right boundary insulated -----tr is not mantained and qr=0
    }                                                                                             // heat convection at right
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
  cout<<'\n';
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
  cout<<'\n';
  cout<<"Heat Flux at right bounadry = "<<hFlux2<<endl;
  }

}
//------------------------------------MAIN PROGRAM------------------------------------------------------------------------------
int main ()
{
  h_inf[0]=0.0;
  h_inf[n-1]=0.0;
  T_inf[0]=0.0;
  T_inf[n-1]=0.0;

  gridGeneration ();
  Gen_Face ();
  Gen ();
  RightBoundary();
  LeftBoundary();
  CoefficientMatrix ();
  sourceVectorCalculation ();	
  matrixInversion ();
 // heatFlux(2.0);            // 1 for left boundary and 2 for right boundary

  cout << "\n";
  cout << " CAQ= ";
  print_array (CAQ);
  cout << "\n"; 
   cout << "\n";
  cout << " S= ";
  print_array (S);
  cout << "\n";

//  print_array(S);
   ofstream graval;
   graval.open("GraphValue.csv",ios::out);
   graval<<"xw,xe,x,dxw,dxe,dx,CAE,CAEO, CAW, CAWO, CAQ, CAQO, CAC, CACO, AT, ATO, CAP, CAPO, S, T"<<endl;

   for(int i=0;i<n;i++)
   {
        graval<<xw[i]<<","<<xe[i]<<","<<x[i]<<","<<dxw[i]<<","<<dxe[i]<<","<<dx[i]<<","<<CAE[i]<<","<<CAEO[i]<<","<< CAW[i]<<","<< CAWO[i]<<","<< CAQ[i]<<","<< CAQO[i]<<","<< CAC[i]<<","<< CACO[i]<<","<< AT[i]<<","<< ATO[i]<<","<< CAP[i]<<","<< CAPO[i]<<","<< S[i]<<","<< T[i]<<endl;
   }
    return 0;;   

  return 0;
}
