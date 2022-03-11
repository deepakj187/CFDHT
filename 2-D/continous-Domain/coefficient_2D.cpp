#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
using namespace std;
/***************************************GLOBAL VARIABLE********************************************************************************/
const int n = 5;

// double T[n+1]={0},S[n*m]={0};
double xw[n] = {0}, xe[n] = {0}, x[n] = {0}, dxw[n] = {0}, dxe[n] = {0}, dx[n] = {0};
const int m = 5;
double ys[m] = {0}, yN[m] = {0}, y[m] = {0}, dys[m] = {0}, dyn[m] = {0}, dy[m] = {0};
// double xy[n][m]={0};
double CAE[n][m] = {0}, CAEO[n][m] = {0}, CAW[n][m] = {0}, CAWO[n][m] = {0}, CAQ[n][m] = {0}, CAQO[n][m] = {0}, CAC[n][m] = {0}, CACO[n][m] = {0}, AT[n][m] = {0}, ATO[n][m] = {0}, CAP[n][m] = {0}, CAPO[n][m] = {0};
double CAS[n][m], CAN[n][m], CASO[n][m], CANO[n][m], S[n][m], T_old[n][m],T[n][m];
double Th_Cond[n][m] = {0.0}, Th_Cond_Face[n + 1][m + 1] = {0.0};

double lx = 0.1; // width of slab            INPUT
double sx = lx / (2.0 * (n - 1));
double ly = 0.18; // width of slab            INPUT
double sy = ly / (2.0 * (m - 1));

double psi = 1.0; // implicit
double mf = 0.0;  // steady state
//******************************************************** temperature *****************
double tl = 0.0;  // left
double tr = 20.0; // right
double tp = 0.0;  // top
double tb = 20.0; // bottom
//+++++++++++++++++++++++++++++++++++++++++++++++ heat flux ++++++++++++++++++++++++++++
double q_l = -10000.0;
double q_r = 0.0;
double q_t = 0.0;
double q_b = 0.0;
//*********************************************** convective ***************************
double h_inf_l = 0.00;
double h_inf_r = 0.00;
double h_inf_t = 20.00;
double h_inf_b = 0.00;

double t_inf_l = 0.00;
double t_inf_r = 0.00;
double t_inf_t = 100.00;
double t_inf_b = 0.00;
//**************************************************************************************Thermo physical properties----------//

double q = 40; // heat generation
double Cp = 1000.0;
double dt = 1.0;
double rho = 2000;

//****************************************************************************************// Boundary conditions variable ***********************//
double l = 1;
double r = 1;
double t = 1;
double b = 1;

double t_l = 1;
double b_l = 1;
double t_r = 1;
double b_r = 1;

/*****************************************************************************************************************************************/

void printArray(double a[])
{

  for (int i = 0; i < m; i++)
    cout << a[i] << " ,";
  // cout<<endl;
}
void print_array(double a[][m])
{
  cout << "\n";
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {

      cout << a[j][i] << " ,";
      // cout<<endl;
    }
    cout << endl;
  }
}
template <typename T, size_t n>
void print_arrayX(T const (&arr)[n])
{
  for (size_t i = 0; i < n; i++)
  {
    std::cout << arr[i] << ',';
  }
}
//-----------------------------------------------------thermoPhysical properties--------------------------------------------------------
void Gen()
{ // defining Area and Th_Cond     INPUT
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      if (i == 0)
      {
        Th_Cond[i][j] = 10.0;
        // q[i][j]=40;
      }
      else if (i == n - 1)
      {
        Th_Cond[i][j] = 10.0;
        // q[i][j]=40;
      }
      else
      {
        Th_Cond[i][j] = 10.0;
        // q[i][j]=40;
      }
    }
  }
}

void Gen_Face()
{ // deefining Area_Face and Th_Cond_Face     INPUT
  for (int i = 0; i <= n; i++)
    for (int j = 0; j < m + 1; j++)
    {
      Th_Cond_Face[i][j] = 10.0;
    }
}

//----------------------------------------------------------GRID GENERATION------------------------------------------------------------------------
void gridGeneration()
{
  xw[0] = 0; // xw[]
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

  xe[0] = sx; // xe[]
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

  x[0] = 0; // x[]
  for (int i = 1; i < n; i++)
  {
    x[i] = x[i - 1] + 2 * sx;
  }

  for (int i = 0; i < n; i++) // dxe dxw dx
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
  print_arrayX(xw);
  cout << "\n";
  cout << "\n";
  cout << " xe= ";
  print_arrayX(xe);
  cout << "\n";
  cout << "\n";
  cout << " x= ";
  print_arrayX(x);
  cout << "\n";
  cout << "\n";
  cout << "dxw= ";
  print_arrayX(dxw);
  cout << "\n";
  cout << "\n";
  cout << "dxe= ";
  print_arrayX(dxe);
  cout << "\n";
  cout << "\n";
  cout << "dx= ";
  print_arrayX(dx);
  cout << "\n";
  ys[0] = 0; // xw[]
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

  yN[0] = sy; // xe[]
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

  y[0] = 0; // x[]
  for (int i = 1; i < m; i++)
  {
    y[i] = y[i - 1] + 2 * sy;
  }

  for (int i = 0; i < m; i++) // dxe dxw dx
  {
    if (i == 0)
    {
      dys[i] = 0;
      dy[i] = sy;
      dyn[i] = 2 * sy;
    }
    else if (i == m - 1)
    {
      dys[i] = 2 * sy;
      dy[i] = sy;
      dyn[i] = 0;
    }
    else
    {
      dys[i] = 2 * sy;
      dyn[i] = 2 * sy;
      dy[i] = 2 * sy;
    }
  }
  cout << "\n";

  cout << " ys= ";
  printArray(ys);
  cout << "\n";
  cout << "\n";
  cout << " yN= ";
  printArray(yN);
  cout << "\n";
  cout << "\n";
  cout << " y= ";
  printArray(y);
  cout << "\n";
  cout << "\n";
  cout << "dys= ";
  printArray(dys);
  cout << "\n";
  cout << "\n";
  cout << "dyn= ";
  printArray(dyn);
  cout << "\n";
  cout << "\n";
  cout << "dy= ";
  printArray(dy);
  cout << "\n";
}

//--------------------------------COEFFICIENT MATRIX------------------------------------------------------------------------------

void CoefficientMatrix()
{
  for (int i = 1; i < n - 1; i++)
  {
    for (int j = 1; j < m - 1; j++)
    {
      CAW[i][j] = (psi * dy[j] * Th_Cond_Face[i][j]) / dxw[i];
      CAE[i][j] = (psi * dy[j] * Th_Cond_Face[i + 1][j]) / dxe[i];
      CAS[i][j] = (psi * dx[i] * Th_Cond_Face[i][j]) / dys[j];
      CAN[i][j] = (psi * dx[i] * Th_Cond_Face[i + 1][j]) / dyn[j];
      CAC[i][j] = 0;
    }
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CAQ[i][j] = (psi * (dx[i] * dy[i] * q));
      AT[i][j] = (mf * (rho * dy[j] * Cp * dx[i]) / dt);
      ATO[i][j] = (mf * (rho * dy[j] * Cp * dx[i]) / dt);
      // S[i][j]=CAQ[i][j] + ATO[i][j]*T[i][j] ;
    }
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      CAP[i][j] = (AT[i][j] + CAW[i][j] + CAE[i][j] + CAC[i][j] + CAS[i][j] + CAN[i][j]);
      CAPO[i][j] = (ATO[i][j] - (CAWO[i][j] + CAEO[i][j] + CASO[i][j] + CANO[i][j])); // for intermediate nodes only---------------------------------
      S[i][j] = CAQ[i][j] + ATO[i][j] * T_old[i][j];
    }
  }
  cout << "\n";
  cout << " CAW= ";
  print_array(CAW);
  cout << "\n";
  cout << "\n";
  cout << " CAE= ";
  print_array(CAE);
  cout << "\n";
  cout << "\n";
  cout << " CAP= ";
  print_array(CAP);
  cout << "\n";
  cout << "\n";
  cout << " CAN= ";
  print_array(CAN);
  cout << "\n";
  cout << "\n";
  cout << " CAS= ";
  print_array(CAS);
  cout << "\n";
  cout << " CAQ= ";
  print_array(CAQ);
  cout << "\n";
  // cout<<"Th_Cond_FACE";
  // for(int i=0;i<=n;i++){
  //   for(int j=0;j<=m;j++){
  //     cout<<Th_Cond_Face[i][j];
  //     cout<<"\n";
  //   }
  // }
}
void BCs_left()
{ // intermediate boundary CV

  for (int j = 1; j < m - 1; j++)
  {
    if (l == 1)
    { // temperature specified
      CAS[0][j] = 0;
      CAN[0][j] = 0;
      CAE[0][j] = 0;
      CAW[0][j] = 0;
      CAP[0][j] = 1;
      S[0][j] = tl;
    }
    if (l == 2)
    { // heat flux specified
      CAS[0][j] = (psi * dx[0] * Th_Cond_Face[0][j]) / dys[j];
      CAN[0][j] = (psi * dx[0] * Th_Cond_Face[0][j + 1]) / dyn[j];
      CAE[0][j] = (psi * dy[j] * Th_Cond_Face[0 + 1][j]) / dxe[0];
      CAW[0][j] = 0;
      CAP[0][j] = CAS[0][j] + CAN[0][j] + CAE[0][j] + CAW[0][j];
      S[0][j] = S[0][j] + dy[j] * q_l * 1.00;
    }
    else
    { // convective specified
      CAS[0][j] = (psi * dx[0] * Th_Cond_Face[0][j]) / dys[j];
      CAN[0][j] = (psi * dx[0] * Th_Cond_Face[0][j + 1]) / dyn[j];
      CAE[0][j] = (psi * dy[j] * Th_Cond_Face[0 + 1][j]) / dxe[0];
      CAW[0][j] = 0;
      CAP[0][j] = CAS[0][j] + CAN[0][j] + CAE[0][j] + CAW[0][j] + dy[j] * 1 * h_inf_l;
      S[0][j] = S[0][j] + dy[j] * h_inf_l * t_inf_l * 1.00;
    }
  }
}
void BCs_right()
{ // intermediate boundary CV

  for (int j = 1; j < m - 1; j++)
  {
    if (r == 1)
    { // temperature specified
      CAS[n - 1][j] = 0;
      CAN[n - 1][j] = 0;
      CAE[n - 1][j] = 0;
      CAW[n - 1][j] = 0;
      CAP[n - 1][j] = 1;
      S[n - 1][j] = tr;
    }
    if (r == 2)
    { // heat flux specified
      CAS[n - 1][j] = (psi * dx[n - 1] * Th_Cond_Face[n - 1][j]) / dys[j];
      CAN[n - 1][j] = (psi * dx[n - 1] * Th_Cond_Face[n - 1][j + 1]) / dyn[j];
      CAE[n - 1][j] = (psi * dy[j] * Th_Cond_Face[n - 1 + 1][j]) / dxe[n - 1];
      CAW[n - 1][j] = 0;
      CAP[n - 1][j] = CAS[n - 1][j] + CAN[n - 1][j] + CAE[n - 1][j] + CAW[n - 1][j];
      S[n - 1][j] = S[n - 1][j] + dy[j] * q_r * 1.00;
    }
    else
    { // convective specified
      CAS[n - 1][j] = (psi * dx[n - 1] * Th_Cond_Face[n - 1][j]) / dys[j];
      CAN[n - 1][j] = (psi * dx[n - 1] * Th_Cond_Face[n - 1][j + 1]) / dyn[j];
      CAE[n - 1][j] = (psi * dy[j] * Th_Cond_Face[n - 1 + 1][j]) / dxe[0];
      CAW[n - 1][j] = 0;
      CAP[n - 1][j] = CAS[n - 1][j] + CAN[n - 1][j] + CAE[n - 1][j] + CAW[n - 1][j] + dy[j] * 1 * h_inf_r;
      S[n - 1][j] = S[n - 1][j] + dy[j] * h_inf_r * t_inf_r * 1.00;
    }
  }
}
void BCs_bottom()
{ // intermediate boundary CV

  for (int i = 1; i < n - 1; i++)
  {
    if (b == 1)
    { // temperature specified
      CAS[i][0] = 0;
      CAN[i][0] = 0;
      CAE[i][0] = 0;
      CAW[i][0] = 0;
      CAP[i][0] = 1;
      S[i][0] = tb;
    }
    if (b == 2)
    { // heat flux specified
      CAS[i][0] = (psi * dx[i] * Th_Cond_Face[i][0]) / dys[0];
      CAN[i][0] = (psi * dx[i] * Th_Cond_Face[i][0 + 1]) / dyn[0];
      CAE[i][0] = (psi * dy[0] * Th_Cond_Face[i + 1][0]) / dxe[i];
      CAW[i][0] = 0;
      CAP[i][0] = CAS[i][0] + CAN[i][0] + CAE[i][0] + CAW[i][0];
      S[i][0] = S[i][0] + dx[i] * q_b * 1.00;
    }
    else
    { // convective specified
      CAS[i][0] = (psi * dx[0] * Th_Cond_Face[i][0]) / dys[0];
      CAN[i][0] = (psi * dx[0] * Th_Cond_Face[i][0 + 1]) / dyn[0];
      CAE[i][0] = (psi * dy[0] * Th_Cond_Face[0 + 1][0]) / dxe[i];
      CAW[i][0] = 0;
      CAP[i][0] = CAS[i][0] + CAN[i][0] + CAE[i][0] + CAW[i][0] + dx[i] * 1 * h_inf_b;
      S[i][0] = S[i][0] + dx[i] * h_inf_b * t_inf_b * 1.00;
    }
  }
}
void BCs_top()
{ // intermediate boundary CV

  for (int i = 1; i < n - 1; i++)
  {
    if (t == 1)
    { // temperature specified
      CAS[i][m - 1] = 0;
      CAN[i][m - 1] = 0;
      CAE[i][m - 1] = 0;
      CAW[i][m - 1] = 0;
      CAP[i][m - 1] = 1;
      S[i][m - 1] = tb;
    }
    if (t == 2)
    { // heat flux specified
      CAS[i][m - 1] = (psi * dx[i] * Th_Cond_Face[i][m - 1]) / dys[m - 1];
      CAN[i][m - 1] = (psi * dx[i] * Th_Cond_Face[i][m - 1 + 1]) / dyn[m - 1];
      CAE[i][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[i + 1][m - 1]) / dxe[i];
      CAW[i][m - 1] = 0;
      CAP[i][m - 1] = CAS[i][m - 1] + CAN[i][m - 1] + CAE[i][m - 1] + CAW[i][m - 1];
      S[i][m - 1] = S[i][m - 1] + dx[i] * q_b * 1.00;
    }
    else
    { // convective specified
      CAS[i][m - 1] = (psi * dx[i] * Th_Cond_Face[i][m - 1]) / dys[m - 1];
      CAN[i][m - 1] = (psi * dx[i] * Th_Cond_Face[i][m - 1 + 1]) / dyn[m - 1];
      CAE[i][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[i + 1][m - 1]) / dxe[i];
      CAW[i][m - 1] = 0;
      CAP[i][m - 1] = CAS[i][m - 1] + CAN[i][m - 1] + CAE[i][m - 1] + CAW[i][m - 1] + dx[i] * 1 * h_inf_b;
      S[i][m - 1] = S[i][m - 1] + dx[i] * h_inf_b * t_inf_b * 1.00;
    }
  }
}
void Corners_tl()
{
  if (t_l = 1) // temp on both boundary
  {
    CAS[0][m - 1] = 0;
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = 0;
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = 1;
    S[0][m - 1] = tl;
  }
  if (t_l = 2) // left=temp-----top=flux
  {
    CAS[0][m - 1] = 0;
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = 0;
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = 1;
    S[0][m - 1] = tl;
  }
  if (t_l = 3) // left=temp-----top=ambient
  {
    CAS[0][m - 1] = 0;
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = 0;
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = 1;
    S[0][m - 1] = tl;
  }
  if (t_l = 4) // left=flux-----top=temp
  {
    CAS[0][m - 1] = 0;
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = 0;
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = 1;
    S[0][m - 1] = tp;
  }
  if (t_l = 5) // left=amb-----top=temp
  {
    CAS[0][m - 1] = 0;
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = 0;
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = 1;
    S[0][m - 1] = tp;
  }
  if (t_l = 6) // left=flux-----top=flux
  {
    CAS[0][m - 1] = (psi * dx[0] * Th_Cond_Face[0][m - 1]) / dys[m - 1];
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[0 + 1][m - 1]) / dxe[0];
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = CAS[0][m - 1] + CAN[0][m - 1] + CAE[0][m - 1] + CAW[0][m - 1];
    S[0][m - 1] = S[0][m - 1] + dxe[0] * q_t * 1.00 + dys[m-1]*q_l*1.00;
  }
   if (t_l = 7) // left=flux-----top=amb
  {
    CAS[0][m - 1] = (psi * dx[0] * Th_Cond_Face[0][m - 1]) / dys[m - 1];
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[0 + 1][m - 1]) / dxe[0];
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = CAS[0][m - 1] + CAN[0][m - 1] + CAE[0][m - 1] + CAW[0][m - 1] +(dxe[0]*h_inf_t);
    S[0][m - 1] = S[0][m - 1] + dxe[0] *h_inf_t*t_inf_t + dys[m-1]*q_l*1.00;
  }
   if (t_l = 8) // left=amb-----top=flux
  {
    CAS[0][m - 1] = (psi * dx[0] * Th_Cond_Face[0][m - 1]) / dys[m - 1];
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[0 + 1][m - 1]) / dxe[0];
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = CAS[0][m - 1] + CAN[0][m - 1] + CAE[0][m - 1] + CAW[0][m - 1] + (dys[m-1]*h_inf_l);
    S[0][m - 1] = S[0][m - 1] + dxe[0] * q_t * 1.00 + dys[m-1]*h_inf_t*t_inf_t;
  }
   if (t_l = 9) // left=amb-----top=amb
  {
    CAS[0][m - 1] = (psi * dx[0] * Th_Cond_Face[0][m - 1]) / dys[m - 1];
    CAN[0][m - 1] = 0;
    CAE[0][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[0 + 1][m - 1]) / dxe[0];
    CAW[0][m - 1] = 0;
    CAP[0][m - 1] = CAS[0][m - 1] + CAN[0][m - 1] + CAE[0][m - 1] + CAW[0][m - 1] + dxe[0] * h_inf_t  + dys[m-1]*h_inf_l;
    S[0][m - 1] = S[0][m - 1] + dxe[0] * h_inf_t*t_inf_t  + dys[m-1]*h_inf_l*t_inf_l;
  }
}
void Corners_bl()
{
  if (b_l = 1) // temp on both boundary
  {
    CAS[0][0] = 0;
    CAN[0][0] = 0;
    CAE[0][0] = 0;
    CAW[0][0] = 0;
    CAP[0][0] = 1;
    S[0][0] = tl;
  }
  if (b_l = 2) // left=temp-----bottom=flux
  {
    CAS[0][0] = 0;
    CAN[0][0] = 0;
    CAE[0][0] = 0;
    CAW[0][0] = 0;
    CAP[0][0] = 1;
    S[0][0] = tl;
  }
  if (b_l = 3) // left=temp-----bottom=ambient
  {
    CAS[0][0] = 0;
    CAN[0][0] = 0;
    CAE[0][0] = 0;
    CAW[0][0] = 0;
    CAP[0][0] = 1;
    S[0][0] = tl;
  }
  if (b_l = 4) // left=flux-----bottom=temp
  {
    CAS[0][0] = 0;
    CAN[0][0] = 0;
    CAE[0][0] = 0;
    CAW[0][0] = 0;
    CAP[0][0] = 1;
    S[0][0] = tb;
  }
  if (b_l = 5) // left=amb-----bottom=temp
  {
    CAS[0][0] = 0;
    CAN[0][0] = 0;
    CAE[0][0] = 0;
    CAW[0][0] = 0;
    CAP[0][0] = 1;
    S[0][0] = tb;
  }
  if (b_l = 6) // left=flux-----bottom=flux
  {
    CAN[0][0] = (psi * dx[0] * Th_Cond_Face[0][0+1]) / dyn[0];
    CAS[0][0] = 0;
    CAE[0][0] = (psi * dy[0] * Th_Cond_Face[0 + 1][0]) / dxe[0];
    CAW[0][0] = 0;
    CAP[0][0] = CAS[0][0] + CAN[0][0] + CAE[0][0] + CAW[0][0];
    S[0][0] = S[0][0] + dxe[0] * q_t * 1.00 + dyn[0]*q_l*1.00;
  }
   if (b_l = 7) // left=flux-----bottom=amb
  {
    CAN[0][0] = (psi * dx[0] * Th_Cond_Face[0][0+1]) / dyn[0];
    CAS[0][0] = 0;
    CAE[0][0] = (psi * dy[0] * Th_Cond_Face[0 + 1][0]) / dxe[0];
    CAW[0][0] = 0;
    CAP[0][0] = CAS[0][0] + CAN[0][0] + CAE[0][0] + CAW[0][0] +(dxe[0]*h_inf_t);
    S[0][0] = S[0][0] + dxe[0] *h_inf_t*t_inf_t + dyn[0]*q_l*1.00;
  }
   if (b_l = 8) // left=amb-----bottom=flux
  {
    CAN[0][0] = (psi * dx[0] * Th_Cond_Face[0][0+1]) / dyn[0];
    CAS[0][0] = 0;
    CAE[0][0] = (psi * dy[0] * Th_Cond_Face[0 + 1][0]) / dxe[0];
    CAW[0][0] = 0;
    CAP[0][0] = CAS[0][0] + CAN[0][0] + CAE[0][0] + CAW[0][0] + (dyn[0]*h_inf_l);
    S[0][0] = S[0][0] + dxe[0] * q_t * 1.00 + dyn[0]*h_inf_t*t_inf_t;
  }
   if (b_l = 9) // left=amb-----bottom=amb
  {
    CAN[0][0] = (psi * dx[0] * Th_Cond_Face[0][0+1]) / dyn[0];
    CAS[0][0] = 0;
    CAE[0][0] = (psi * dy[0] * Th_Cond_Face[0 + 1][0]) / dxe[0];
    CAW[0][0] = 0;
    CAP[0][0] = CAS[0][0] + CAN[0][0] + CAE[0][0] + CAW[0][0] + dxe[0] * h_inf_t  + dyn[0]*h_inf_l;
    S[0][0] = S[0][0] + dxe[0] * h_inf_t*t_inf_t  + dyn[0]*h_inf_l*t_inf_l;
  }
}
void Corners_tr()
{
  if (t_r = 1) // temp on both boundary
  {
    CAS[n-1][m - 1] = 0;
    CAN[n-1][m - 1] = 0;
    CAE[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = 1;
    S[n-1][m - 1] = tr;
  }
  if (t_r = 2) // right=temp-----top=flux
  {
    CAS[n-1][m - 1] = 0;
    CAN[n-1][m - 1] = 0;
    CAE[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = 1;
    S[n-1][m - 1] = tr;
  }
  if (t_r = 3) // right=temp-----top=ambient
  {
    CAS[n-1][m - 1] = 0;
    CAN[n-1][m - 1] = 0;
    CAE[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = 1;
    S[n-1][m - 1] = tr;
  }
  if (t_r = 4) // right=flux-----top=temp
  {
    CAS[n-1][m - 1] = 0;
    CAN[n-1][m - 1] = 0;
    CAE[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = 1;
    S[n-1][m - 1] = tp;
  }
  if (t_r = 5) // right=amb-----top=temp
  {
    CAS[n-1][m - 1] = 0;
    CAN[n-1][m - 1] = 0;
    CAE[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = 1;
    S[n-1][m - 1] = tp;
  }
  if (t_r = 6) // right=flux-----top=flux
  {
    CAS[n-1][m - 1] = (psi * dx[n-1] * Th_Cond_Face[n-1][m - 1]) / dys[m - 1];
    CAN[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[n-1][m - 1]) / dxw[n-1];
    CAE[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = CAS[n-1][m - 1] + CAN[n-1][m - 1] + CAE[n-1][m - 1] + CAW[n-1][m - 1];
    S[n-1][m - 1] = S[n-1][m - 1] + dxw[n-1] * q_t * 1.00 + dys[m-1]*q_r*1.00;
  }
   if (t_r = 7) // right=flux-----top=amb
  {
    CAS[n-1][m - 1] = (psi * dx[n-1] * Th_Cond_Face[n-1][m - 1]) / dys[m - 1];
    CAN[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[n-1][m - 1]) / dxw[n-1];
    CAE[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = CAS[n-1][m - 1] + CAN[n-1][m - 1] + CAE[n-1][m - 1] + CAW[n-1][m - 1] +(dxw[n-1]*h_inf_t);
    S[n-1][m - 1] = S[n-1][m - 1] + dxw[n-1] *h_inf_t*t_inf_t + dys[m-1]*q_l*1.00;
  }
   if (t_r = 8) // right=amb-----top=flux
  {
    CAS[n-1][m - 1] = (psi * dx[n-1] * Th_Cond_Face[n-1][m - 1]) / dys[m - 1];
    CAN[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[n-1][m - 1]) / dxw[n-1];
    CAE[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = CAS[n-1][m - 1] + CAN[n-1][m - 1] + CAE[n-1][m - 1] + CAW[n-1][m - 1] + (dys[m-1]*h_inf_l);
    S[n-1][m - 1] = S[n-1][m - 1] + dxw[n-1] * q_t * 1.00 + dys[m-1]*h_inf_t*t_inf_t;
  }
   if (t_r = 9) // right=amb-----top=amb
  {
    CAS[n-1][m - 1] = (psi * dx[n-1] * Th_Cond_Face[n-1][m - 1]) / dys[m - 1];
    CAN[n-1][m - 1] = 0;
    CAW[n-1][m - 1] = (psi * dy[m - 1] * Th_Cond_Face[n-1][m - 1]) / dxw[n-1];
    CAE[n-1][m - 1] = 0;
    CAP[n-1][m - 1] = CAS[n-1][m - 1] + CAN[n-1][m - 1] + CAE[n-1][m - 1] + CAW[n-1][m - 1] + dxw[n-1] * h_inf_t  + dys[m-1]*h_inf_l;
    S[n-1][m - 1] = S[n-1][m - 1] + dxw[n-1] * h_inf_t*t_inf_t  + dys[m-1]*h_inf_l*t_inf_l;
  }
}
void Corners_br()
{
  if (b_r = 1) // temp on both boundary
  {
    CAS[n-1][0] = 0;
    CAN[n-1][0] = 0;
    CAE[n-1][0] = 0;
    CAW[n-1][0] = 0;
    CAP[n-1][0] = 1;
    S[n-1][0] = tr;
  }
  if (b_r = 2) // right=temp-----bottom=flux
  {
    CAS[n-1][0] = 0;
    CAN[n-1][0] = 0;
    CAE[n-1][0] = 0;
    CAW[n-1][0] = 0;
    CAP[n-1][0] = 1;
    S[n-1][0] = tr;
  }
  if (b_r = 3) // right=temp-----bottom=ambient
  {
    CAS[n-1][0] = 0;
    CAN[n-1][0] = 0;
    CAE[n-1][0] = 0;
    CAW[n-1][0] = 0;
    CAP[n-1][0] = 1;
    S[n-1][0] = tr;
  }
  if (b_r = 4) // right=flux-----bottom=temp
  {
    CAS[n-1][0] = 0;
    CAN[n-1][0] = 0;
    CAE[n-1][0] = 0;
    CAW[n-1][0] = 0;
    CAP[n-1][0] = 1;
    S[n-1][0] = tb;
  }
  if (b_r = 5) // right=amb-----bottom=temp
  {
    CAS[n-1][0] = 0;
    CAN[n-1][0] = 0;
    CAE[n-1][0] = 0;
    CAW[n-1][0] = 0;
    CAP[n-1][0] = 1;
    S[n-1][0] = tb;
  }
  if (b_r = 6) // right=flux-----bottom=flux
  {
    CAN[n-1][0] = (psi * dx[n-1] * Th_Cond_Face[n-1][0+1]) / dyn[0];
    CAS[n-1][0] = 0;
    CAW[n-1][0] = (psi * dy[n-1] * Th_Cond_Face[n-1][0]) / dxw[n-1];
    CAE[n-1][0] = 0;
    CAP[n-1][0] = CAS[n-1][0] + CAN[n-1][0] + CAE[n-1][0] + CAW[n-1][0];
    S[n-1][0] = S[n-1][0] + dxw[n-1] * q_b * 1.00 + dyn[0]*q_r*1.00;
  }
   if (b_r = 7) // right=flux-----bottom=amb
  {
    CAN[n-1][0] = (psi * dx[n-1] * Th_Cond_Face[n-1][0+1]) / dyn[0];
    CAS[n-1][0] = 0;
    CAW[n-1][0] = (psi * dy[n-1] * Th_Cond_Face[n-1][0]) / dxw[n-1];
    CAE[n-1][0] = 0;
    CAP[n-1][0] = CAS[n-1][0] + CAN[n-1][0] + CAE[n-1][0] + CAW[n-1][0] +(dxw[n-1]*h_inf_b);
    S[n-1][0] = S[n-1][0] + dxw[n-1] *h_inf_b*t_inf_b + dyn[0]*q_r*1.00;
  }
   if (b_r = 8) // right=amb-----bottom=flux
  {
    CAN[n-1][0] = (psi * dx[n-1] * Th_Cond_Face[n-1][0+1]) / dyn[0];
    CAS[n-1][0] = 0;
    CAW[n-1][0] = (psi * dy[n-1] * Th_Cond_Face[n-1][0]) / dxw[n-1];
    CAE[n-1][0] = 0;
    CAP[n-1][0] = CAS[n-1][0] + CAN[n-1][0] + CAE[n-1][0] + CAW[n-1][0] + (dyn[0]*h_inf_r);
    S[n-1][0] = S[n-1][0] + dxw[n-1] * q_b * 1.00 + dyn[0]*h_inf_r*t_inf_r;
  }
   if (b_r = 9) // right=amb-----bottom=amb
  {
    CAN[n-1][0] = (psi * dx[n-1] * Th_Cond_Face[n-1][0+1]) / dyn[0];
    CAS[n-1][0] = 0;
    CAW[n-1][0] = (psi * dy[n-1] * Th_Cond_Face[n-1][0]) / dxw[n-1];
    CAE[n-1][0] = 0;
    CAP[n-1][0] = CAS[n-1][0] + CAN[n-1][0] + CAE[n-1][0] + CAW[n-1][0] + dxw[n-1] * h_inf_b  + dyn[0]*h_inf_r;
    S[n-1][0] = S[n-1][0] + dxw[n-1] * h_inf_b*t_inf_b  + dyn[0]*h_inf_r*t_inf_r;
  }
}
void ADIscheme(){

} 

int main()
{

  Gen();
  Gen_Face();
  gridGeneration();
  CoefficientMatrix();
  BCs_left();
  BCs_right();
  BCs_top();
  BCs_bottom();
  Corners_tl();
  Corners_tr();
  Corners_bl();
  Corners_br();
  ADIscheme();

}