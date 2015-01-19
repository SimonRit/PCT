#include "pctfillholl_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>

#include "SmallHoleFiller.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkChangeInformationImageFilter.h>

int main(int argc, char * argv[])
{
  GGO(pctfillholl, args_info);

  typedef float OutputPixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  // Reader
  typedef itk::ImageFileReader< OutputImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( args_info.input_arg );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( reader->Update() );

  SmallHoleFiller< OutputImageType > filler;
  filler.SetImage( reader->GetOutput() );
  filler.SetHolePixel(0.);
  filler.Fill();

  typedef itk::ChangeInformationImageFilter< OutputImageType > CIIType;
  CIIType::Pointer cii = CIIType::New();
  cii->SetInput(filler.GetOutput());
  cii->ChangeOriginOn();
  cii->ChangeDirectionOn();
  cii->ChangeSpacingOn();
  cii->SetOutputDirection( reader->GetOutput()->GetDirection() );
  cii->SetOutputOrigin(    reader->GetOutput()->GetOrigin() );
  cii->SetOutputSpacing(   reader->GetOutput()->GetSpacing() );

  // Write
  typedef itk::ImageFileWriter<  OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( cii->GetOutput() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );

  return EXIT_SUCCESS;
}
