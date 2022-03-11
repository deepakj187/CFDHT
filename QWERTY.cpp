#include<iostream>
#include<algorithm>
#include<cmath>
#include<fstream>
using namespace std;
/*********GLOBAL VARIABLE*********/
const int n=11;                                               //// pls change the array size if it is not 5      INPUT
double A[n]={0},B[n]={0},T[n+1]={0},S[n]={0},r=0.002; 
double Xw[n]={0}, Xe[n]={0}, X[n]={0}, dxw[n]={0}, dxe[n]={0}, dx[n]={0};
double  q=250000.0, Pp[n]={0.0}, h_inf[n]={0.0}, rho[n]={0.0}, Cp=0.0, dt=1.0, T_inf[n]={0.0};                     //INPUT
double CAE[n]={0}, CAEO[n]={0}, CAW[n]={0}, CAWO[n]={0}, CAQ[n]={0}, CAQO[n]={0}, CAC[n]={0}, CACO[n]={0}, AT[n]={0}, ATO[n]={0}, CAP[n]={0}, CAPO[n]={0};
double Area[n]={0.0},Area_Face[n+1]={0.0},Th_Cond[n]={0.0},Th_Cond_Face[n+1]={0.0};

/**********************/

template<typename T, size_t n>
void print_array(T const(& arr)[n])
{
    for (size_t i = 0; i < n; i++) {
        std::cout << arr[i] << ',';
    }
}
void Gen(){                                                         //defining Area and Th_Cond     INPUT
    for(int i=0;i<n;i++)
    {
     Area[i]=4*atan(1)*(r*r);
     Th_Cond[i]=5.0;
     Pp[i]=2*4*atan(1)*r;
     h_inf[i]=12.5;
     T_inf[i]=20;
    }
}

void Gen_Face(){                                                     //deefining Area_Face and Th_Cond_Face     INPUT
    for(int i=0;i<=n;i++)
    { 
        Area_Face[i]=4*atan(1)*(r*r);
        Th_Cond_Face[i]=5.0;
    }
}

void gridGeneration(double s, int n)                                            //Grid Generation
{
   Xw[0]=0;                        // xw[]
    for (int i=1;i<n;i++)
    {
        if (i==1)
        {
            Xw[i]=Xw[i-1] + s;
        }
        else
        {
            Xw[i]=Xw[i-1] + 2*s;
        }
    }
    
    Xe[0]=s;                        // xe[]
    for (int i=1;i<n;i++)
    {
        if (i==n-1)
        {
            Xe[i]=Xe[i-1] + s;
        }
        else
        {
            Xe[i]=Xe[i-1] + 2*s;
        }
    }
    
    X[0]=0;                         // x[]
    for (int i=1;i<n;i++)
    {
        X[i]=X[i-1] + 2*s;
        
    }
    
    for(int i=0;i<n;i++)        // dxe dxw dx
    {
        if (i==0 )
        {
            dxw[i]=0;
            dx[i]=s;
            dxe[i]=2*s;
        }
        else if (i==n-1)
        {
            dxw[i]=2*s;
            dx[i]=s;
            dxe[i]=0;
        }
        else
        {
            dxw[i]=2*s;
            dxe[i]=2*s;
            dx[i]=2*s;
        }
    }
    cout<<"\n";
    
    cout<<" xw= ";
    print_array(Xw);
    cout<<"\n";
    cout<<"\n";
    cout<<" xe= ";
    print_array(Xe);
    cout<<"\n";
    cout<<"\n";
    cout<<" x= ";
    print_array(X);
    cout<<"\n";
    cout<<"\n";
    cout<<"dxw= ";
    print_array(dxw);
    cout<<"\n";
    cout<<"\n";
    cout<<"dxe= ";
    print_array(dxe);
    cout<<"\n";
    cout<<"\n";
    cout<<"dx= ";
    print_array(dx);
    cout<<"\n";
    
}


void CoefficientMatrix(int n ,double psi,double mf)                             //Coefficient Matrix
{
    CAE[0]=0;                                   //BCs       INPUT
    CAW[0]=0;
    CAE[n-1]=0;                                   //BCs       INPUT
    CAW[n-1]=0;
    CAP[0]=1;                                   //BCs       INPUT
    CAP[n-1]=1;                                   //BCs       INPUT
    
    for(int i=1;i<n-1;i++)
    {
        CAE[i]=(psi*Area_Face[i+1]*Th_Cond_Face[i+1])/dxe[i];
        CAEO[i]=((1-psi)*Area_Face[i+1]*Th_Cond_Face[i+1])/dxe[i];
        CAW[i]=(psi*Area_Face[i]*Th_Cond_Face[i])/dxw[i];
        CAWO[i]=((1-psi)*Area_Face[i]*Th_Cond_Face[i])/dxw[i];
        CAQ[i]=(psi*(Area[i]*q*dx[i]));
        CAQO[i]=(1-psi)*(Area[i]*q*dx[i]);
        CAC[i]=(psi*(Pp[i]*h_inf[i]*dx[i]));
        CACO[i]=(1-psi)*(Pp[i]*h_inf[i]*dx[i]);
        AT[i]=(mf*(rho[i]*Area[i]*Cp*dx[i])/dt);
        ATO[i]=(mf*(rho[i]*Area[i]*Cp*dx[i])/dt);
    }
    for(int i=1;i<n-1;i++)
    {
        CAP[i]=(AT[i]+CAW[i]+CAE[i]+CAC[i]);
        CAPO[i]=(ATO[i]-(CAWO[i]+CAEO[i]+CACO[i]));

    }
    cout<<"\n";
    cout<<" CAW= ";
    print_array(CAW);
    cout<<"\n";
    cout<<"\n";
    cout<<" CAE= ";
    print_array(CAE);
    cout<<"\n";
    cout<<"\n";
    cout<<" CAP= ";
    print_array(CAP);
    cout<<"\n";
   
}

void sourceVectorCalculation(double t1,double tm){                       //sourceVectorCalculation(double t0,double tn)
    S[0]=t1;     //left boundary
    S[n-1]=tm;     // right boundary
    for (int i=1;i<n-1;i++)
        S[i]= ((CAPO[i]*T[i]) + (CAWO[i]*T[i-1]) + (CAEO[i]*T[i+1]) + CAQO[i] + CAQ[i] + CAC[i]*T_inf[i] + CACO[i]*T_inf[i]);
}

void matrixInversion(){                                 //matrixInversion()
// forward substitution
    for (int i=0;i<n;i++){
        if(i==0){
        A[i]=CAE[i]/CAP[i];
        B[i]=S[i]/CAP[i];
        }
        else
        {
            A[i]=CAE[i]/(CAP[i]-CAW[i]*A[i-1]);
            B[i]=(S[i]+CAW[i]*B[i-1])/(CAP[i]-CAW[i]*A[i-1]);
        }
    }
    cout<<"\n";
    cout<<"[A]= ";
    print_array(A);
    cout<<"\n";
    cout<<"\n";
    cout<<"[B]= ";
    print_array(B);
    cout<<"\n";
//  back substitution

    for (int i=n-1;i>=0;i--){                           
        
        T[i]=A[i]*T[i+1] + B[i];
    }
    cout<<"\n";
    cout<<"the desired solution are [T]=";
    for (int i=0;i<n;i++){
        cout<<T[i]<<",";
    }
    cout<<"\n";
}

int main()
{
    double l=0.040;                                  // width of slab            INPUT
    int n=11;                                        // no of nodes
    double s = l / (2.0 * (n - 1));
    
    gridGeneration(s,n);
    Gen_Face();
    Gen();
    
    CoefficientMatrix(n,1.00,0.00);                     //n,shi,MF              INPUT
    sourceVectorCalculation(60,80);             //t1 && tm//            INPUT
    matrixInversion();
   ofstream graval;
   graval.open("GraphValue.csv",ios::out);
   graval<<"Xw,Xe,X,dxw,dxe,dx,CAE,CAEO, CAW, CAWO, CAQ, CAQO, CAC, CACO, AT, ATO, CAP, CAPO, S, T"<<endl;

   for(int i=0;i<n;i++)
   {
        graval<<Xw[i]<<","<<Xe[i]<<","<<X[i]<<","<<dxw[i]<<","<<dxe[i]<<","<<dx[i]<<","<<CAE[i]<<","<<CAEO[i]<<","<< CAW[i]<<","<< CAWO[i]<<","<< CAQ[i]<<","<< CAQO[i]<<","<< CAC[i]<<","<< CACO[i]<<","<< AT[i]<<","<< ATO[i]<<","<< CAP[i]<<","<< CAPO[i]<<","<< S[i]<<","<< T[i]<<endl;
   }
    return 0;
}