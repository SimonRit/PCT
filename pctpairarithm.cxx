#include "pctpairarithm_ggo.h"

#include <rtkMacro.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

int main(int argc, char * argv[])
{
  GGO(pctpairarithm, args_info);

  // Read
  typedef itk::Vector<float, 3> VectorType;
  typedef itk::Image<VectorType,2> PairsImageType;
  typedef itk::ImageFileReader< PairsImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(args_info.input_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( reader->Update() );

  if(args_info.vector_given != 3)
    {
    std::cerr << "ERROR: you must pass a 3D vector in --vector" << std::endl;
    return EXIT_FAILURE;
    }

  // Multiply
  itk::ImageRegionIteratorWithIndex<PairsImageType> it(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    VectorType v = it.Get();
    if(it.GetIndex()[0]!=4)
      {
      v[0] *= args_info.vector_arg[0];
      v[1] *= args_info.vector_arg[1];
      v[2] *= args_info.vector_arg[2];
      }
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
