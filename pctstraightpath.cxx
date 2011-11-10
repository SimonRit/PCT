#include "pctstraightpath_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>
#include <itkRayQuadricIntersectionFunction.h>
#include <itkConstantImageSource.h>
#include <itkFDKBackProjectionImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDivideImageFilter.h>
#include <itkImageIterator.h>

#include "pctBetheBlochFunctor.h"

#define PAIRS_IN_RAM 1000000

int main(int argc, char * argv[])
{
  GGO(pctstraightpath, args_info);

  // Create a stack of empty projection images
  typedef float OutputPixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::ConstantImageSource< OutputImageType > ConstantImageSourceType;
  ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctstraightpath>(constantImageSource, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( constantImageSource->Update() );

  // Create corresponding stack to count events
  typedef itk::Image< unsigned int, Dimension > CountImageType;
  typedef itk::ConstantImageSource< CountImageType > CountImageSourceType;
  CountImageSourceType::Pointer countImageSource = CountImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<CountImageSourceType, args_info_pctstraightpath>(countImageSource, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( countImageSource->Update() );

  // Create pipeline for normalization and output
  typedef itk::DivideImageFilter<OutputImageType, CountImageType, OutputImageType> DivideImageFilterType;
  DivideImageFilterType::Pointer div = DivideImageFilterType::New();
  div->SetInput1(constantImageSource->GetOutput());
  div->SetInput2(countImageSource->GetOutput());

  typedef itk::ImageFileWriter<  OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( div->GetOutput() );

  typedef itk::ImageFileWriter<  CountImageType > CountWriterType;
  CountWriterType::Pointer cwriter = CountWriterType::New();
  cwriter->SetFileName( args_info.count_arg );
  cwriter->SetInput( countImageSource->GetOutput() );

  // Read pairs
  typedef itk::Vector<float, Dimension> VectorType;
  typedef itk::Image<VectorType,2> PairsImageType;
  typedef itk::ImageFileReader< PairsImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( args_info.input_arg );
  reader->UpdateOutputInformation();

  // Image information constants
  OutputImageType::SizeType imgSize  = constantImageSource->GetOutput()->GetBufferedRegion().GetSize();
  OutputImageType::PointType imgOrigin = constantImageSource->GetOutput()->GetOrigin();
  OutputImageType::SpacingType imgSpacing = constantImageSource->GetOutput()->GetSpacing();
  float *imgData = constantImageSource->GetOutput()->GetBufferPointer();
  unsigned int *imgCountData = countImageSource->GetOutput()->GetBufferPointer();
  unsigned long npixelsPerSlice = imgSize[0] * imgSize[1];
  itk::Vector<float, Dimension> imgSpacingInv;
  OutputImageType::PointType originInVox;
  for(unsigned int i=0; i<Dimension; i++)
    {
    imgSpacingInv[i] = 1./imgSpacing[i];
    originInVox[i] = -imgOrigin[i]*imgSpacingInv[i];
    }

  size_t nprotons = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];

  pct::Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double> convFunc;
  PairsImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();

  unsigned int nregions = nprotons/PAIRS_IN_RAM+1;


  // Corrections
  typedef itk::Vector<double,3> DVectorType;
  DVectorType source;
  source.Fill(0.);
  source[2] = args_info.source_arg;

  typedef itk::RayQuadricIntersectionFunction<double,3> RQIType;
  RQIType::Pointer inSphere = RQIType::New();
  inSphere->SetA(1.);
  inSphere->SetB(1.);
  inSphere->SetC(1.);
  inSphere->SetI(-2*args_info.source_arg);
  inSphere->SetJ(args_info.source_arg*args_info.source_arg - 855*855);

  RQIType::Pointer outSphere = RQIType::New();
  outSphere->SetA(1.);
  outSphere->SetB(1.);
  outSphere->SetC(1.);
  outSphere->SetI(-2*args_info.source_arg);
  outSphere->SetJ(args_info.source_arg*args_info.source_arg - 1145*1145);

  // Create magnitude lut, need one proton pair to know plane positions
  region.SetIndex(1, 0);
  region.SetSize(1, 1);
  reader->GetOutput()->SetRequestedRegion(region);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( reader->Update() );

  std::vector<double> zmag(imgSize[2]);
  const double sourcePosInVox = (args_info.source_arg-imgOrigin[2]) * imgSpacingInv[2];
  const double zPlaneOutInMM = (*(reader->GetOutput()->GetBufferPointer()+1))[2];
  const double zPlaneOutInVox = (zPlaneOutInMM-imgOrigin[2]) * imgSpacingInv[2];
  for(unsigned int i=0; i<imgSize[2]; i++)
    zmag[i] = (zPlaneOutInVox-sourcePosInVox)/(i-sourcePosInVox);

  for(unsigned int r=0; r<nregions; r++)
  {
    // Read r-th set of pairs
    region.SetIndex(1, r*PAIRS_IN_RAM);
    region.SetSize(1, vnl_math_min(PAIRS_IN_RAM, int(nprotons-region.GetIndex(1))));
    reader->GetOutput()->SetRequestedRegion(region);
    TRY_AND_EXIT_ON_ITK_EXCEPTION( reader->Update() );

    // Process pairs
    itk::ImageRegionIterator<PairsImageType> it(reader->GetOutput(), region);
    for(size_t p=region.GetIndex(1); p<region.GetIndex(1)+region.GetSize(1); p++)
    {
      if(p%PAIRS_IN_RAM==0)
        {
        std::cout << '\r' 
                  << p << " pairs of protons processed ("
                  << 100*p/nprotons << "%), writing current results"
                  << std::flush;
        TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );
        if(args_info.count_given)
          TRY_AND_EXIT_ON_ITK_EXCEPTION( cwriter->Update() );
        }
      DVectorType pIn = it.Get();
      ++it;
      DVectorType pOut = it.Get();
      ++it;
      DVectorType dIn = it.Get();
      ++it;
      DVectorType dOut = it.Get();
      ++it;
      const double value = convFunc.GetValue(it.Get()[1], it.Get()[0]);
      ++it;

      // Move straight to entrance and exit shapes
      inSphere->SetRayOrigin(pIn);
      outSphere->SetRayOrigin(pOut);
      if(!inSphere->Evaluate(dIn) || !outSphere->Evaluate(dOut))
        //Weird angle, probably a nuclear interaction
        continue;

      DVectorType pSIn  = pIn  + dIn *inSphere->GetFarthestDistance();
      if(pSIn[2]<args_info.source_arg)
        pSIn  = pIn  + dIn * inSphere->GetNearestDistance();
      DVectorType pSOut = pOut + dOut*outSphere->GetFarthestDistance();
      if(pSOut[2]<args_info.source_arg)
        pSOut  = pOut  + dOut * outSphere->GetNearestDistance();

      // Convert everything to voxel coordinates
      for(unsigned int i=0; i<Dimension; i++)
        {
        pSIn[i]  = (pSIn[i]  - imgOrigin[i]) * imgSpacingInv[i];
        pSOut[i] = (pSOut[i] - imgOrigin[i]) * imgSpacingInv[i];
        pIn[i]   = (pIn[i]   - imgOrigin[i]) * imgSpacingInv[i];
        pOut[i]  = (pOut[i]  - imgOrigin[i]) * imgSpacingInv[i];
        dIn[i]   = dIn[i]  * imgSpacingInv[i];
        dOut[i]  = dOut[i] * imgSpacingInv[i];
        }

      // Normalize direction with respet to z
      dIn[0] /= dIn[2];
      dIn[1] /= dIn[2];
      //dIn[2] = 1.; SR: implicit in the following
      dOut[0] /= dOut[2];
      dOut[1] /= dOut[2];
      //dOut[2] = 1.; SR: implicit in the following

      // Parameters of the 3rd order polynomial. The function goes from
      // z=0 to z=lastSliceZ
      double ax = pSIn[0];
      double ay = pSIn[1];
      double x  = pSOut[0];
      double y  = pSOut[1];
      const double bx = dIn[0];
      const double by = dIn[1];
      const double xd = dOut[0];
      const double yd = dOut[1];
      const double lastSliceZ = (pSOut[2]-pSIn[2]);
      const double invzsq = 1./(lastSliceZ*lastSliceZ);
      const double cx = invzsq * ( 3*x - lastSliceZ*xd - 3*ax - 2*bx*lastSliceZ );
      const double cy = invzsq * ( 3*y - lastSliceZ*yd - 3*ay - 2*by*lastSliceZ );
      const double inv3zsq = invzsq/3.;
      const double dx = inv3zsq * ( xd - bx - 2*cx*lastSliceZ );
      const double dy = inv3zsq * ( yd - by - 2*cy*lastSliceZ );

      for(unsigned int k=0; k<imgSize[2]; k+=1)
        {
        double xx, yy;
        const double dk = k; //SR make sure to convert it once only
        if(dk<=pSIn[2])
          {
          const double z = (dk-pIn[2]);
          xx = pIn[0]+z*dIn[0];
          yy = pIn[1]+z*dIn[1];
          }
        else if(dk>=pSOut[2])
          {
          const double z = (dk-pSOut[2]);
          xx = pSOut[0]+z*dOut[0];
          yy = pSOut[1]+z*dOut[1];
          }
        else
          {
          const double z = (dk-pSIn[2]);
          xx = ax+z*(bx+z*(cx+z*dx));
          yy = ay+z*(by+z*(cy+z*dy));
          }

        xx = (xx-originInVox[0])*zmag[k]+originInVox[0];
        const int i = int(xx+0.5);
        if(i<0 || i>=(int)imgSize[0])
          continue;

        yy = (yy-originInVox[1])*zmag[k]+originInVox[1];
        const int j = int(yy+0.5);
        if(j<0 || j>=(int)imgSize[1])
          continue;

        unsigned long idx = i+j*imgSize[0]+k*npixelsPerSlice;
        imgData[ idx ] += value;
        imgCountData[ idx ]++;
        }
      constantImageSource->GetOutput()->Modified();
      countImageSource->GetOutput()->Modified();
    }
  }
  std::cout << "Finalizing binning and writing..." << std::endl;
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );
  if(args_info.count_given)
    TRY_AND_EXIT_ON_ITK_EXCEPTION( cwriter->Update() );

  return EXIT_SUCCESS;
}
