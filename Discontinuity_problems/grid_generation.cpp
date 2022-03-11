#include <iostream>
#include <stdio.h>
#include<iomanip>

#include<algorithm>
#include<cmath>
#include<fstream>
using namespace std;

using namespace std;
//-------------------------------------------------------------------INPUTS------------------------------------------
int const d=3;
double l[d+1]={10,5,5,20};
int M[d+1]={5,5,5,5};
double s[d+1]={0};
const int n=15;
//--------------------------------------------------------------------------------------------------------------------
template < typename T, size_t n > void
print_array (T const (&arr)[n]){
  for (size_t i = 0; i < n; i++)
    {
      std::cout << arr[i] << ',';
    }
}



double xw[n]={0}, xe[n]={0}, x[n]={0}, dxw[n]={0}, dxe[n]={0}, dx[n]={0},li[d-1];
//-----------------------------------------------------GRID GENERATION------------------------------------------------
void gridgeneration()
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
    
    //-----------------------------------------------------------RESULTS------------------------------------------------------------
    cout<<setw(10)<<"xw"<<setw(10)<<"xe"<<setw(10)<<"x"<<endl;
    for (int i=0;i<n;i++)
    {cout<<"i="<<i+1<<setw(8)<<xw[i]<<setw(10)<<xe[i]<<setw(10)<<x[i]<<endl;
    }
    cout<<setw(10)<<"dxw"<<setw(10)<<"dxe"<<setw(10)<<"dx"<<endl;
    for (int i=0;i<n;i++)
    {cout<<"i="<<i+1<<setw(8)<<dxw[i]<<setw(10)<<dxe[i]<<setw(10)<<dx[i]<<endl;
    }//-------------------------------------------------------------------------------------------------------------------------------
}

int main()
{
    gridgeneration();
    cout<<"\n";
    print_array(s);
}