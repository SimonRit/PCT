#include "pctpaircuts_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>
#include <itkConstantImageSource.h>

#include "itkProtonPairsToDistanceDrivenProjection.h"

#include <itkImageFileWriter.h>
#include <itkRegularExpressionSeriesFileNames.h>
#include <itkTimeProbe.h>

#define PAIRS_IN_RAM 1000000

int main(int argc, char * argv[])
{
  GGO(pctpaircuts, args_info);

  typedef float ProjectionPixelType;
  typedef itk::Image< ProjectionPixelType, 2 > ProjectionImageType;

  // Create a stack of empty projection images and associated proton count
  typedef itk::ConstantImageSource< ProjectionImageType > ConstantImageSourceType;
  ConstantImageSourceType::Pointer sumEnergy   = ConstantImageSourceType::New();
  ConstantImageSourceType::Pointer sumEnergySq = ConstantImageSourceType::New();
  ConstantImageSourceType::Pointer sumAngle    = ConstantImageSourceType::New();
  ConstantImageSourceType::Pointer sumAngleSq  = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumEnergy, args_info);
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumEnergySq, args_info);
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumAngle, args_info);
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumAngleSq, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumEnergy->Update() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumEnergySq->Update() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumAngle->Update() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumAngleSq->Update() );

  typedef itk::Image< unsigned int, 2 > CountImageType;
  typedef itk::ConstantImageSource< CountImageType > CountImageSourceType;
  CountImageSourceType::Pointer counts = CountImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<CountImageSourceType, args_info_pctpaircuts>(counts, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( counts->Update() );

  // Read pairs
  typedef itk::Vector<float, 3> VectorType;
  typedef itk::Image<VectorType,2> PairsImageType;
  typedef itk::ImageFileReader< PairsImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( args_info.input_arg );
  reader->UpdateOutputInformation();
  size_t nprotons = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
  PairsImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
  unsigned int nregions = nprotons/PAIRS_IN_RAM+1;

  // Image information constants
  const ProjectionImageType::SizeType imgSize  = sumEnergy->GetOutput()->GetBufferedRegion().GetSize();
  const ProjectionImageType::PointType imgOrigin = sumEnergy->GetOutput()->GetOrigin();
  const ProjectionImageType::SpacingType imgSpacing = sumEnergy->GetOutput()->GetSpacing();
  itk::Vector<float, 2> imgSpacingInv;
  for(unsigned int i=0; i<2; i++)
    imgSpacingInv[i] = 1./imgSpacing[i];

  // Data pointers
  float *pSumEnergy     = sumEnergy->GetOutput()->GetBufferPointer();
  float *pSumEnergySq   = sumEnergySq->GetOutput()->GetBufferPointer();
  float *pSumAngle      = sumAngle->GetOutput()->GetBufferPointer();
  float *pSumAngleSq    = sumAngleSq->GetOutput()->GetBufferPointer();
  unsigned int *pCounts = counts->GetOutput()->GetBufferPointer();

  std::cout << "Compute cuts..." << std::endl;
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
      VectorType pIn = it.Get();
      ++it;
      VectorType pOut = it.Get();
      ++it;
      VectorType dIn = it.Get();
      ++it;
      VectorType dOut = it.Get();
      ++it;
      VectorType data = it.Get();
      ++it;

      static double mag = (args_info.source_arg - pOut[2]) / (args_info.source_arg - pIn[2]);

      const double xx = (pIn[0]*mag-imgOrigin[0]) * imgSpacingInv[0];
      const int i = int(xx+0.5);
      if(i<0 || i>=(int)imgSize[0])
        continue;

      const double yy = (pIn[1]*mag-imgOrigin[1]) * imgSpacingInv[1];
      const int j = int(yy+0.5);
      if(j<0 || j>=(int)imgSize[1])
        continue;

      const unsigned long idx = i+j*imgSize[0];
      const double angle = vcl_acos(dIn * dOut);
      const double energy = data[1];
      pSumEnergy  [idx] += energy;
      pSumEnergySq[idx] += energy*energy;
      pSumAngle   [idx] += angle;
      pSumAngleSq [idx] += angle*angle;
      pCounts     [idx]++;
    }
  }

  std::cout << "Finalize cuts..." << std::endl;

  // Now compute the cuts and the average
  size_t npixels = sumEnergy->GetOutput()->GetBufferedRegion().GetNumberOfPixels();
  for(unsigned int idx=0; idx<npixels; idx++)
    {
    if(pCounts[idx])
      {
      pSumEnergy  [idx] /= pCounts[idx];
      pSumAngle   [idx] /= pCounts[idx];
      pSumEnergySq[idx] /= pCounts[idx];
      pSumAngleSq [idx] /= pCounts[idx];
      pSumEnergySq[idx] = sqrt( pSumEnergySq[idx]-pSumEnergy[idx]*pSumEnergy[idx] );
      pSumAngleSq[idx]  = sqrt( pSumAngleSq [idx]-pSumAngle [idx]*pSumAngle [idx] );
      }
    }
  for(unsigned int idx=0; idx<npixels; idx++)
    {
    pSumEnergySq[idx] *= 3;
    pSumAngleSq[idx]  *= 3;
    }

  std::cout << "Select pairs..." << std::endl;

  // And select the pairs
  std::vector<VectorType> pairs;
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
      VectorType pIn = it.Get();
      ++it;
      VectorType pOut = it.Get();
      ++it;
      VectorType dIn = it.Get();
      ++it;
      VectorType dOut = it.Get();
      ++it;
      VectorType data = it.Get();
      ++it;

      static double mag = (args_info.source_arg - pOut[2]) / (args_info.source_arg - pIn[2]);

      const double xx = (pIn[0]*mag-imgOrigin[0]) * imgSpacingInv[0];
      const int i = int(xx+0.5);
      if(i<0 || i>=(int)imgSize[0])
        continue;

      const double yy = (pIn[1]*mag-imgOrigin[1]) * imgSpacingInv[1];
      const int j = int(yy+0.5);
      if(j<0 || j>=(int)imgSize[1])
        continue;

      const unsigned long idx = i+j*imgSize[0];
      const double angle = vcl_acos(dIn * dOut);
      if( vcl_abs(angle  -pSumAngle [idx])<pSumAngleSq [idx] &&
          vcl_abs(data[1]-pSumEnergy[idx])<pSumEnergySq[idx] )
        {
        pairs.push_back(pIn);
        pairs.push_back(pOut);
        pairs.push_back(dIn);
        pairs.push_back(dOut);
        pairs.push_back(data);
        }
    }
  }

  std::cout << "Write pairs..." << std::endl;

  itk::ImageRegion<2> pairsRegion;
  itk::ImageRegion<2>::SizeType size;
  size[0] = 5;
  size[1] = pairs.size()/5;
  pairsRegion.SetSize(size);
  PairsImageType::Pointer img = PairsImageType::New();
  img->SetRegions(pairsRegion);
  img->Allocate();

  itk::ImageRegionIterator<PairsImageType> it(img, pairsRegion);
  for(size_t i=0; i<pairs.size(); ++i, ++it)
    it.Set( pairs[i] );

  // Write
  typedef itk::ImageFileWriter<PairsImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( img );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );

  // Optional outputs
  typedef itk::ImageFileWriter<ProjectionImageType> PWriterType;
  if(args_info.mangle_given)
    {
    PWriterType::Pointer w = PWriterType::New();
    w->SetInput(sumAngle->GetOutput());
    w->SetFileName(args_info.mangle_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }
  if(args_info.sangle_given)
    {
    PWriterType::Pointer w = PWriterType::New();
    w->SetInput(sumAngleSq->GetOutput());
    w->SetFileName(args_info.sangle_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }
  if(args_info.menergy_given)
    {
    PWriterType::Pointer w = PWriterType::New();
    w->SetInput(sumEnergy->GetOutput());
    w->SetFileName(args_info.menergy_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }  typedef itk::ImageFileWriter<ProjectionImageType> PWriterType;
  if(args_info.senergy_given)
    {
    PWriterType::Pointer w = PWriterType::New();
    w->SetInput(sumEnergySq->GetOutput());
    w->SetFileName(args_info.senergy_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }
  if(args_info.count_given)
    {
    itk::ImageFileWriter<CountImageType>::Pointer w;
    w = itk::ImageFileWriter<CountImageType>::New();
    w->SetInput(counts->GetOutput());
    w->SetFileName(args_info.count_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }
  return EXIT_SUCCESS;
}
