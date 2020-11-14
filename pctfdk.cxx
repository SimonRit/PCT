#include "pctfdk_ggo.h"
#include "rtkGgoFunctions.h"

#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkProjectionsReader.h"
#include "pctDDParkerShortScanImageFilter.h"
#include "pctFDKDDConeBeamReconstructionFilter.h"

#include <itkRegularExpressionSeriesFileNames.h>
#include <itkImageFileWriter.h>

int main(int argc, char * argv[])
{
  GGO(pctfdk, args_info);

  typedef float OutputPixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads( std::min<double>(8, itk::MultiThreaderBase::GetGlobalMaximumNumberOfThreads() ) );

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

  // Projections reader
  typedef itk::Image< OutputPixelType, Dimension+1 > ProjectionImageType;
  typedef rtk::ProjectionsReader< ProjectionImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames( names->GetFileNames() );
  if(args_info.wpc_given)
    {
    std::vector<double> coeffs;
    coeffs.assign(args_info.wpc_arg, args_info.wpc_arg+args_info.wpc_given);
    reader->SetWaterPrecorrectionCoefficients(coeffs);
    }

  // Geometry
  if(args_info.verbose_flag)
    std::cout << "Reading geometry information from "
              << args_info.geometry_arg
              << "..."
              << std::endl;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(args_info.geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( geometryReader->GenerateOutputInformation() )

  // Short scan image filter
  typedef pct::DDParkerShortScanImageFilter< ProjectionImageType > PSSFType;
  PSSFType::Pointer pssf = PSSFType::New();
  pssf->SetInput( reader->GetOutput() );
  pssf->SetGeometry( geometryReader->GetOutputObject() );
  pssf->InPlaceOff();

  // Create reconstructed image
  typedef rtk::ConstantImageSource< OutputImageType > ConstantImageSourceType;
  ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctfdk>(constantImageSource, args_info);

  // FDK reconstruction filtering
  typedef pct::FDKDDConeBeamReconstructionFilter< OutputImageType > FDKCPUType;
  FDKCPUType::Pointer feldkamp = FDKCPUType::New();
  feldkamp->SetInput( 0, constantImageSource->GetOutput() );
  feldkamp->SetProjectionStack( pssf->GetOutput() );
  feldkamp->SetGeometry( geometryReader->GetOutputObject() );
  feldkamp->GetRampFilter()->SetTruncationCorrection(args_info.pad_arg);
  feldkamp->GetRampFilter()->SetHannCutFrequency(args_info.hann_arg);
  feldkamp->GetRampFilter()->SetHannCutFrequencyY(args_info.hannY_arg);                                                       \

  // Write
  typedef itk::ImageFileWriter<  OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( feldkamp->GetOutput() );

  if(args_info.verbose_flag)
    std::cout << "Reconstructing and writing... " << std::flush;
  itk::TimeProbe writerProbe;

  writerProbe.Start();
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );
  writerProbe.Stop();

  if(args_info.verbose_flag)
    {
    std::cout << "It took " << writerProbe.GetMean() << ' ' << writerProbe.GetUnit() << std::endl;
    feldkamp->PrintTiming(std::cout);
    }

  return EXIT_SUCCESS;
}
