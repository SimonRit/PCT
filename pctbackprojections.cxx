#include "pctbackprojections_ggo.h"
#include "rtkGgoFunctions.h"

#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "pctFDKDDBackProjectionImageFilter.h"
#include "rtkJosephBackProjectionImageFilter.h"

#include <itkRegularExpressionSeriesFileNames.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>
#include <itkImageSeriesReader.h>

int main(int argc, char * argv[])
{
  GGO(pctbackprojections, args_info);

  typedef float OutputPixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image< OutputPixelType, Dimension+1 > ProjectionImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  // Geometry
  if(args_info.verbose_flag)
    std::cout << "Reading geometry information from "
              << args_info.geometry_arg
              << "..."
              << std::flush;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(args_info.geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( geometryReader->GenerateOutputInformation() )
  if(args_info.verbose_flag)
    std::cout << " done." << std::endl;

  // Create an empty volume
  typedef rtk::ConstantImageSource< OutputImageType > ConstantImageSourceType;
  ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctbackprojections>(constantImageSource, args_info);

  // Generate file names
  itk::RegularExpressionSeriesFileNames::Pointer names = itk::RegularExpressionSeriesFileNames::New();
  names->SetDirectory(args_info.path_arg);
  names->SetNumericSort(false);
  names->SetRegularExpression(args_info.regexp_arg);
  names->SetSubMatch(0);

  if(args_info.verbose_flag)
    std::cout << "Reading "
              << names->GetFileNames().size()
              << " projection file(s)..."
              << std::flush;

  // Projections reader
  typedef itk::ImageSeriesReader< ProjectionImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames( names->GetFileNames() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( reader->Update() );
  if(args_info.verbose_flag)
    std::cout << " done." << std::endl;

  // Create back projection image filter
  if(args_info.verbose_flag)
    std::cout << "Backprojecting volume..." << std::flush;
  typedef pct::FDKDDBackProjectionImageFilter<OutputImageType, OutputImageType> BPType;
  BPType::Pointer bp = BPType::New();
  bp->SetInput( constantImageSource->GetOutput() );
  bp->SetProjectionStack( reader->GetOutput() );
  bp->SetGeometry( geometryReader->GetOutputObject() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( bp->Update() )
  if(args_info.verbose_flag)
    std::cout << " done ." << std::endl;

  // Write
  if(args_info.verbose_flag)
    std::cout << "Writing... " << std::flush;
  typedef itk::ImageFileWriter<  OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( bp->GetOutput() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );
  if(args_info.verbose_flag)
    std::cout << " done" << std::endl;

  return EXIT_SUCCESS;
}
