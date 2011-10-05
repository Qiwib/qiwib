#include "FDInterval.cpp"

class Gaussian : public FDInterval::RealFunction {
public:
  double x0, a;

  Gaussian(double x0, double a) : x0(x0), a(a) {}

  double operator()(double x) const {
    return exp(-a*(x-x0)*(x-x0));
  }
};

class GaussianCurvature : public FDInterval::RealFunction {
public:
  double x0, a;

  GaussianCurvature(double x0, double a) : x0(x0), a(a) {}

  double operator()(double x) const {
    double r2 = (x-x0)*(x-x0);
    return sqrt(M_E/(32.0*a))*fabs((2.0*a*r2-1.0))*2.0*a*exp(-a*r2);
  }
};


int main()
{
  vector<double> grid(11);
  for(size_t i=0;i<=10;i++) grid[i] = (i-5.0)/3.0;
  FDInterval space(grid);
  double exponent = 100.0;
  vector<double> centers(1,0.0);
  Gaussian blip(centers[0],exponent);
  GaussianCurvature blip_curvature(centers[0],exponent);

  space.refine_around_points(centers,5.0,1.0);

  cout << "coarsedensity = " << FDInterval::FDFunction(space,blip) << ";\n";

  space.refine_to_function(blip_curvature,0,0.005);
  FDInterval::FDFunction density(space, blip);

  cout << "density   = " << density << ";\n";
  cout << "curvature = " << FDInterval::FDFunction(space, blip_curvature) << ";\n";

  FDInterval::FDFunction potential(space.solve_poisson(blip));

  cout << "potential = " << potential << ";\n";
  return 0;
}
