#include<iostream>
#include<algorithm>
#include<cmath>
#include<fstream>
using namespace std;
/**************************GLOBAL VARIABLE************************/
const int n=5;
const int m=10; 
//double T[n+1]={0},S[n*m]={0}; 
double xw[n]={0}, xe[n]={0}, x[n]={0}, dxw[n]={0}, dxe[n]={0}, dx[n]={0};      
double ys[m]={0}, yn[m]={0}, y[m]={0}, dys[m]={0}, dyN[m]={0}, dy[m]={0};

double lx= 0.1;		// width of slab            INPUT
double sx = lx / (2.0 * (n - 1));
double ly = 0.18;		// width of slab            INPUT
double sy = ly / (2.0 * (m - 1));

double psi=1.0;                            //implicit 
double mf=0.0;                             // steady state 


/******************************************************************/

void printArray(double a[]){
	for(int i=0;i<m;i++)
      cout<< a[i]<<" ,";
  	//cout<<endl;
}
template < typename T, size_t n > void
print_arrayX (T const (&arr)[n]){
  for (size_t i = 0; i < n; i++)
    {
      std::cout << arr[i] << ',';
    }
}


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
  print_arrayX (xw);
  cout << "\n";
  cout << "\n";
  cout << " xe= ";
  print_arrayX (xe);
  cout << "\n";
  cout << "\n";
  cout << " x= ";
  print_arrayX (x);
  cout << "\n";
  cout << "\n";
  cout << "dxw= ";
  print_arrayX (dxw);
  cout << "\n";
  cout << "\n";
  cout << "dxe= ";
  print_arrayX (dxe);
  cout << "\n";
  cout << "\n";
  cout << "dx= ";
  print_arrayX (dx);
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

  yn[0] = sy;			// xe[]
  for (int i = 1; i < m; i++)
    {
      if (i == m - 1)
	{
	  yn[i] = yn[i - 1] + sy;
	}
      else
	{
	  yn[i] = yn[i - 1] + 2 * sy;
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
  printArray (ys);
  cout << "\n";
  cout << "\n";
  cout << " yn= ";
  printArray (yn);
  cout << "\n";
  cout << "\n";
  cout << " y= ";
  printArray (y);
  cout << "\n";
  cout << "\n";
  cout << "dys= ";
  printArray (dys);
  cout << "\n";
  cout << "\n";
  cout << "dyN= ";
  printArray (dyN);
  cout << "\n";
  cout << "\n";
  cout << "dy= ";
  printArray (dy);
  cout << "\n";

}
int main(){
    gridGeneration();
}
