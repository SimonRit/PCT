#include "pctzengbackprojections_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>

#include "pctZengBackProjectionImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int main(int argc, char * argv[])
{
  GGO(pctzengbackprojections, args_info);

  typedef float InputPixelType;
  const unsigned int Dimension = 4;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;

  itk::ImageFileReader<InputImageType>::Pointer reader;
  reader = itk::ImageFileReader<InputImageType>::New();
  reader->SetFileName(args_info.input_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->Update())

  typedef pct::ZengBackProjectionImageFilter<InputImageType> ZengFilterType;
  ZengFilterType::Pointer zeng;
  zeng = ZengFilterType::New();
  zeng->SetInput( reader->GetOutput() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION(zeng->Update())

  typedef itk::ImageFileWriter< ZengFilterType::OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(zeng->GetOutput(0));
  writer->SetFileName(args_info.outputc_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->UpdateLargestPossibleRegion())
  writer->SetInput(zeng->GetOutput(1));
  writer->SetFileName(args_info.outputs_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->UpdateLargestPossibleRegion())


  return EXIT_SUCCESS;
}
