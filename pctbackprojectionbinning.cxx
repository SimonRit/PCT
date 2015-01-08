#include "pctbackprojectionbinning_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>
#include <rtkConstantImageSource.h>
#include <rtkThreeDCircularProjectionGeometryXMLFile.h>

#include "pctProtonPairsToBackProjection.h"
#include "SmallHoleFiller.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRegularExpressionSeriesFileNames.h>
#include <itkTimeProbe.h>
#include <itkChangeInformationImageFilter.h>

int main(int argc, char * argv[])
{
  GGO(pctbackprojectionbinning, args_info);

  typedef float OutputPixelType;
  const unsigned int Dimension = 4;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef pct::ProtonPairsToBackProjection<OutputImageType, OutputImageType> ProjectionFilter;

  // Generate file names
  itk::RegularExpressionSeriesFileNames::Pointer names = itk::RegularExpressionSeriesFileNames::New();
  names->SetDirectory(args_info.path_arg);
  names->SetNumericSort(false);
  names->SetRegularExpression(args_info.regexp_arg);
  names->SetSubMatch(0);

  if(args_info.verbose_flag)
    std::cout << "Regular expression matches "
              << names->GetFileNames().size()
              << " file(s)..."
              << std::endl;

  // Read or create bp images
  OutputImageType::Pointer inBp;
  ProjectionFilter::CountImagePointer inCount;
  if(args_info.bpVal_given)
    {
    if(args_info.verbose_flag)
      std::cout << "Reading input files "
                << args_info.bpVal_arg
                << " and "
                << args_info.bpCount_arg
                << "..."
                << std::endl;

    itk::ImageFileReader<OutputImageType>::Pointer readBPVal = itk::ImageFileReader<OutputImageType>::New();
    readBPVal->SetFileName(args_info.bpVal_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(readBPVal->Update());
    inBp = readBPVal->GetOutput();

    itk::ImageFileReader<ProjectionFilter::CountImageType>::Pointer readBPCount = itk::ImageFileReader<ProjectionFilter::CountImageType>::New();
    readBPCount->SetFileName(args_info.bpCount_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(readBPCount->Update());
    inCount = readBPCount->GetOutput();
    }
  else
    {
    if(args_info.verbose_flag)
      std::cout << "Starting from empty images..."
                << std::endl;
    rtk::ConstantImageSource< OutputImageType >::Pointer constantImageSource = rtk::ConstantImageSource< OutputImageType >::New();
    rtk::SetConstantImageSourceFromGgo<rtk::ConstantImageSource< OutputImageType >, args_info_pctbackprojectionbinning>(constantImageSource, args_info);
    TRY_AND_EXIT_ON_ITK_EXCEPTION( constantImageSource->Update() );
    inBp = constantImageSource->GetOutput();

    rtk::ConstantImageSource< ProjectionFilter::CountImageType >::Pointer countSource = rtk::ConstantImageSource< ProjectionFilter::CountImageType >::New();
    rtk::SetConstantImageSourceFromGgo<rtk::ConstantImageSource< ProjectionFilter::CountImageType >, args_info_pctbackprojectionbinning>(countSource, args_info);
    TRY_AND_EXIT_ON_ITK_EXCEPTION( countSource->Update() );
    inCount = countSource->GetOutput();
    }

  // Projection filter
  ProjectionFilter::Pointer projection = ProjectionFilter::New();
  projection->SetInput( inBp );
  projection->SetCount( inCount );
  projection->SetProtonPairsFileNames( names->GetFileNames() );
  projection->SetMostLikelyPathType( args_info.mlptype_arg );
  projection->SetIonizationPotential( args_info.ionpot_arg * CLHEP::eV );

  // Geometry
  if(args_info.verbose_flag)
    std::cout << "Reading geometry information from "
              << args_info.geometryIn_arg
              << " and "
              << args_info.geometryOut_arg
              << "..."
              << std::endl;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReaderIn;
  geometryReaderIn = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReaderIn->SetFilename(args_info.geometryIn_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( geometryReaderIn->GenerateOutputInformation() )
  projection->SetGeometryIn( geometryReaderIn->GetOutputObject() );
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReaderOut;
  geometryReaderOut = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReaderOut->SetFilename(args_info.geometryOut_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( geometryReaderOut->GenerateOutputInformation() )
  projection->SetGeometryOut( geometryReaderOut->GetOutputObject() );

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

  // Write
  typedef itk::ImageFileWriter<  OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( cii->GetOutput() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );

  if(args_info.count_given)
    {
    // Write
    typedef itk::ImageFileWriter< ProjectionFilter::CountImageType > CountWriterType;
    CountWriterType::Pointer cwriter = CountWriterType::New();
    cwriter->SetFileName( args_info.count_arg );
    cwriter->SetInput( projection->GetCount() );
    TRY_AND_EXIT_ON_ITK_EXCEPTION( cwriter->Update() )
    }

  return EXIT_SUCCESS;
}
