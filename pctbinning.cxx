#include "pctbinning_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>
#include <rtkConstantImageSource.h>

#include "pctProtonPairsToDistanceDrivenProjection.h"
#include "SmallHoleFiller.h"

#include <itkImageFileWriter.h>
#include <itkRegularExpressionSeriesFileNames.h>
#include <itkTimeProbe.h>
#include <itkChangeInformationImageFilter.h>

int main(int argc, char * argv[])
{
  GGO(pctbinning, args_info);

  if(args_info.elosswepl_given && args_info.output_given)
    {
    std::cerr << "Only --output or --elosswepl should be provided" << std::endl;
    return EXIT_FAILURE;
    }

  typedef float OutputPixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

#if ( ( ITK_VERSION_MAJOR > 4 ) )
  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads( std::min<double>(8., itk::MultiThreaderBase::GetGlobalMaximumNumberOfThreads() ) );
#else
  itk::MultiThreader::SetGlobalMaximumNumberOfThreads( std::min<double>(8., itk::MultiThreader::GetGlobalMaximumNumberOfThreads() ) );
#endif

  // Create a stack of empty projection images
  typedef rtk::ConstantImageSource< OutputImageType > ConstantImageSourceType;
  ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctbinning>(constantImageSource, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( constantImageSource->Update() );

  // Projection filter
  typedef pct::ProtonPairsToDistanceDrivenProjection<OutputImageType, OutputImageType> ProjectionFilter;
  ProjectionFilter::Pointer projection = ProjectionFilter::New();
  projection->SetInput( constantImageSource->GetOutput() );
  projection->SetProtonPairsFileName( args_info.input_arg );
  projection->SetSourceDistance( args_info.source_arg );
  projection->SetMostLikelyPathType( args_info.mlptype_arg );
  projection->SetIonizationPotential( args_info.ionpot_arg * CLHEP::eV );
  projection->SetRobust( args_info.robust_flag );
  projection->SetComputeScattering( args_info.scatwepl_given );

  if(args_info.quadricIn_given)
    {
    //quadric = object surface
    ProjectionFilter::RQIType::Pointer qIn = ProjectionFilter::RQIType::New();
    qIn->SetA(args_info.quadricIn_arg[0]);
    qIn->SetB(args_info.quadricIn_arg[1]);
    qIn->SetC(args_info.quadricIn_arg[2]);
    qIn->SetD(args_info.quadricIn_arg[3]);
    qIn->SetE(args_info.quadricIn_arg[4]);
    qIn->SetF(args_info.quadricIn_arg[5]);
    qIn->SetG(args_info.quadricIn_arg[6]);
    qIn->SetH(args_info.quadricIn_arg[7]);
    qIn->SetI(args_info.quadricIn_arg[8]);
    qIn->SetJ(args_info.quadricIn_arg[9]);
    projection->SetQuadricIn(qIn);
    }
  if(args_info.quadricOut_given)
    {
    ProjectionFilter::RQIType::Pointer qOut = ProjectionFilter::RQIType::New();
    qOut->SetA(args_info.quadricOut_arg[0]);
    qOut->SetB(args_info.quadricOut_arg[1]);
    qOut->SetC(args_info.quadricOut_arg[2]);
    qOut->SetD(args_info.quadricOut_arg[3]);
    qOut->SetE(args_info.quadricOut_arg[4]);
    qOut->SetF(args_info.quadricOut_arg[5]);
    qOut->SetG(args_info.quadricOut_arg[6]);
    qOut->SetH(args_info.quadricOut_arg[7]);
    qOut->SetI(args_info.quadricOut_arg[8]);
    qOut->SetJ(args_info.quadricOut_arg[9]);
    projection->SetQuadricOut(qOut);
    }

  TRY_AND_EXIT_ON_ITK_EXCEPTION( projection->Update() );

  SmallHoleFiller< OutputImageType > filler;
  if(args_info.fill_flag)
    {
    filler.SetImage( projection->GetOutput() );
    filler.SetHolePixel(0.);
    filler.Fill();
    }

  typedef itk::ChangeInformationImageFilter< OutputImageType > CIIType;
  CIIType::Pointer cii = CIIType::New();
  if(args_info.fill_flag)
    cii->SetInput(filler.GetOutput());
  else
    cii->SetInput(projection->GetOutput());
  cii->ChangeOriginOn();
  cii->ChangeDirectionOn();
  cii->ChangeSpacingOn();
  cii->SetOutputDirection( projection->GetOutput()->GetDirection() );
  cii->SetOutputOrigin(    projection->GetOutput()->GetOrigin() );
  cii->SetOutputSpacing(   projection->GetOutput()->GetSpacing() );

  if(args_info.elosswepl_given || args_info.output_given)
    {
    // Write
    typedef itk::ImageFileWriter<  OutputImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    if(args_info.elosswepl_given)
      writer->SetFileName( args_info.elosswepl_arg );
    else
      writer->SetFileName( args_info.output_arg );
    writer->SetInput( cii->GetOutput() );
    TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );
    }

  if(args_info.count_given)
    {
    // Write
    typedef itk::ImageFileWriter< ProjectionFilter::CountImageType > CountWriterType;
    CountWriterType::Pointer cwriter = CountWriterType::New();
    cwriter->SetFileName( args_info.count_arg );
    cwriter->SetInput( projection->GetCount() );
    TRY_AND_EXIT_ON_ITK_EXCEPTION( cwriter->Update() )
    }

  if(args_info.scatwepl_given)
    {
    // Write
    typedef itk::ImageFileWriter< ProjectionFilter::AngleImageType > AngleWriterType;
    AngleWriterType::Pointer swriter = AngleWriterType::New();
    swriter->SetFileName( args_info.scatwepl_arg );
    swriter->SetInput( projection->GetAngle() );
    TRY_AND_EXIT_ON_ITK_EXCEPTION( swriter->Update() )
    }

  return EXIT_SUCCESS;
}
