#include <iostream>
using std::cout;
#include <vector>
#include <cmath>	// using sqrt

class Field
{
  private:
    int b_; 
    int nx_;
    int ny_;
    std::vector<double> data_; 
  public:
    Field(int nx, int ny) : 
      b_(nx/10), nx_(nx), ny_(ny), data_  ((nx+1) * (ny+1) ) 
    {
      for (int i=0;i<(nx+1)*(ny+1);i++) data_[i]=0.0;
    }

    double  operator() (int x, int y) const 
    {
      return data_[x + (nx_+1)*y];
    }

    double& operator() (int x, int y) 
    {
      return data_[x + (nx_+1)*y];
    }

    bool on_boundaryA(int x, int y) const 
    {
      if (x>=2*b_ && x<=6*b_ && 
          y>=7*b_ && y<=9*b_) return true;
      return false; 
    }
 
    bool on_boundaryB(int x, int y) const 
    {
      if (x==8*b_ && y>=1*b_ && y<=5*b_) return true; 
      if (y==5*b_ && x>=6*b_ && y<=8*b_) return true; 
      return false; 
    }

// Central difference for dphi/dy
    double dy(double xx, double yy) const 
    {
      int x = nx_*xx;
      int y = ny_*yy;
      return 0.5*double(nx_)*(operator()(x,y+1)-operator()(x,y-1)); 
    }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
};


double sor_update(Field& phi, double omega)
{
  int nx=phi.nx(); 
  int ny=phi.ny(); 
  double d2=0.0;
  for (int x=1;x<nx;x++)
    for (int y=1;y<ny;y++)
    {
      if (phi.on_boundaryA(x,y)) phi(x,y) = 1.0; 
      else if (phi.on_boundaryB(x,y)) phi(x,y) = -1.0; 
      else 
      {
        double phi_old = phi(x,y);
        phi(x,y) = (1.0 - omega) * phi(x,y) + 
          + 0.25 * omega * 
          ( phi(x+1,y  )+ phi(x-1,y  )+ phi(x  ,y+1)+ phi(x  ,y-1) ); 
        double d = phi(x,y) - phi_old; 
        d2 += d*d;
      }
    }
  return sqrt(d2/double(nx*ny));
}

// Fixed number of iterations
void sor(Field& phi, double omega, int niter)
{ 
  for (int iter=0;iter<niter;iter++) 
    sor_update(phi,omega);
}

// Usage: ./q1 <grid width> <number of iterations>
int main(int argc, char* argv[])
{
  if (argc!=3) return 1;
  double omega=1.99;

// Grid resolution - n=100 is good enough for gnuplot
  int n=atoi(argv[1]);
// Fix number of SOR iterations 
  int niter=atoi(argv[2]);

  Field phi(n,n);

  for (int y=0;y<n+1;y++)
    phi(n,y) = y/double(n);
  for (int x=0;x<n+1;x++)
    phi(x,n) = x/double(n);
  
  sor(phi,omega,niter); 
  cout << phi.dy(0.3,0.5) << "\n";

  return 0;
}
