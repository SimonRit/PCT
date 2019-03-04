#include "pctmostlikelypath_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>
#include <rtkQuadricShape.h>

#include "pctSchulteMLPFunction.h"
#include "pctThirdOrderPolynomialMLPFunction.h"
#include "pctPolynomialMLPFunction.h"

#include <itkImageFileWriter.h>

int main(int argc, char * argv[])
{
  GGO(pctmostlikelypath, args_info);

  typedef itk::Vector<double,3> VectorType;
  const unsigned int Dimension = 1;
  typedef itk::Image< VectorType, Dimension > OutputImageType;

  // Create outputs
  OutputImageType::RegionType region;
  region.SetSize(0,args_info.dimension_arg);
  OutputImageType::Pointer trajectory = OutputImageType::New();
  trajectory->SetRegions(region);
  trajectory->Allocate();
  region.SetSize(0,2);
  OutputImageType::Pointer intersections = OutputImageType::New();
  intersections->SetRegions(region);
  intersections->Allocate();

  typedef rtk::QuadricShape RQIType;
  RQIType::Pointer quadricIn = RQIType::New();
  if(args_info.quadricIn_given)
    {
    quadricIn->SetA(args_info.quadricIn_arg[0]);
    quadricIn->SetB(args_info.quadricIn_arg[1]);
    quadricIn->SetC(args_info.quadricIn_arg[2]);
    quadricIn->SetD(args_info.quadricIn_arg[3]);
    quadricIn->SetE(args_info.quadricIn_arg[4]);
    quadricIn->SetF(args_info.quadricIn_arg[5]);
    quadricIn->SetG(args_info.quadricIn_arg[6]);
    quadricIn->SetH(args_info.quadricIn_arg[7]);
    quadricIn->SetI(args_info.quadricIn_arg[8]);
    quadricIn->SetJ(args_info.quadricIn_arg[9]);
    }

  RQIType::Pointer quadricOut = RQIType::New();
  if(args_info.quadricOut_given)
    {
    quadricOut->SetA(args_info.quadricOut_arg[0]);
    quadricOut->SetB(args_info.quadricOut_arg[1]);
    quadricOut->SetC(args_info.quadricOut_arg[2]);
    quadricOut->SetD(args_info.quadricOut_arg[3]);
    quadricOut->SetE(args_info.quadricOut_arg[4]);
    quadricOut->SetF(args_info.quadricOut_arg[5]);
    quadricOut->SetG(args_info.quadricOut_arg[6]);
    quadricOut->SetH(args_info.quadricOut_arg[7]);
    quadricOut->SetI(args_info.quadricOut_arg[8]);
    quadricOut->SetJ(args_info.quadricOut_arg[9]);
    }

  VectorType pIn(args_info.posIn_arg);
  VectorType pOut(args_info.posOut_arg);
  VectorType dIn(args_info.dirIn_arg);
  VectorType dOut(args_info.dirOut_arg);

  // Move straight to entrance and exit shapes
  VectorType pSIn  = pIn;
  VectorType pSOut = pOut;
  double nearDistIn, nearDistOut, farDistIn, farDistOut;
  if(args_info.quadricIn_given &&
     args_info.quadricOut_given &&
     quadricIn->IsIntersectedByRay(pIn,dIn,nearDistIn,farDistIn) &&
     quadricOut->IsIntersectedByRay(pOut,dOut,nearDistOut,farDistOut))
    {
    pSIn  = pIn  + dIn  * nearDistIn;
    if(pSIn[2]<pIn[2]  || pSIn[2]>pOut[2])
      pSIn  = pIn  + dIn  * farDistIn;
    pSOut = pOut + dOut * nearDistOut;
    if(pSOut[2]<pIn[2] || pSOut[2]>pOut[2])
      pSOut = pOut + dOut * farDistOut;
    }
  OutputImageType::IndexType index;
  index[0] = 0;
  intersections->SetPixel(index,pSIn);
  index[0] = 1;
  intersections->SetPixel(index,pSOut);

  // Normalize direction with respect to z
  dIn[0] /= dIn[2];
  dIn[1] /= dIn[2];
  //dIn[2] = 1.; SR: implicit in the following
  dOut[0] /= dOut[2];
  dOut[1] /= dOut[2];
  //dOut[2] = 1.; SR: implicit in the following

  //pct::ThirdOrderPolynomialMLPFunction<double>::Pointer mlp;
  pct::MostLikelyPathFunction<double>::Pointer mlp;
  if(args_info.type_arg==std::string("schulte"))
    mlp = pct::SchulteMLPFunction::New();
  else if(args_info.type_arg==std::string("polynomial"))
    mlp = pct::ThirdOrderPolynomialMLPFunction<double>::New();
  else if(args_info.type_arg==std::string("krah"))
    {
      pct::PolynomialMLPFunction::Pointer mlp_poly;
      mlp_poly = pct::PolynomialMLPFunction::New();
      // pct::PolynomialMLPFunction::Pointer polynomial_mlp = dynamic_cast<pct::PolynomialMLPFunction*>(mlp.GetPointer());
      mlp_poly->SetPolynomialDegree(args_info.mlppolydeg_arg);
      mlp = mlp_poly;
    }
  else
    {
    std::cerr << "Unknown mlp type: " << args_info.type_arg << std::endl;
    exit(1);
    }
  mlp->Init(pSIn, pSOut, dIn, dOut);

  for(int k=0; k<args_info.dimension_arg; k++)
    {
    double x, y;
    const double u = args_info.origin_arg+k*args_info.spacing_arg;
    if(u<=pSIn[2])
      {
      const double z = (u-pIn[2]);
      x = pIn[0]+z*dIn[0];
      y = pIn[1]+z*dIn[1];
      }
    else if(u>=pSOut[2])
      {
      const double z = (u-pSOut[2]);
      x = pSOut[0]+z*dOut[0];
      y = pSOut[1]+z*dOut[1];
      }
    else
      {
      mlp->Evaluate(u, x, y);
      }
    VectorType point;
    point[0] = x;
    point[1] = y;
    point[2] = u;

    index[0] = k;
    trajectory->SetPixel(index, point);
    }


  // Write
  typedef itk::ImageFileWriter<  OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( args_info.trajectory_arg );
  writer->SetInput( trajectory );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );

  writer->SetFileName( args_info.intersections_arg );
  writer->SetInput( intersections );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );

  return EXIT_SUCCESS;
}
