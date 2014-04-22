#include "pctswapcoordinates_ggo.h"

#include <rtkMacro.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

int main(int argc, char * argv[])
{
  GGO(pctswapcoordinates, args_info);

  // Read
  typedef itk::Vector<float, 3> VectorType;
  typedef itk::Image<VectorType,2> PairsImageType;
  typedef itk::ImageFileReader< PairsImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(args_info.input_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( reader->Update() );

  // Swap
  itk::ImageRegionIteratorWithIndex<PairsImageType> it(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    VectorType v = it.Get();
    if(it.GetIndex()[0]!=4)
      std::swap(v[(args_info.fixed_arg+1)%3], v[(args_info.fixed_arg+2)%3]);
    it.Set(v);
    }

  // Write
  typedef itk::ImageFileWriter< PairsImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(args_info.output_arg);
  writer->SetInput( reader->GetOutput() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );

  return EXIT_SUCCESS;
}
