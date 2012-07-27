#include "itkImageFileWriter.h" // For intermediate debugging output
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkNumericTraits.h"
#include "rtkMacro.h"

template <typename TImage>
SmallHoleFiller<TImage>::SmallHoleFiller()
{
  this->HolePixel = itk::NumericTraits< typename TImage::PixelType >::Zero;
  this->Image = NULL;
  this->Output = NULL;
}

template <typename TImage>
void SmallHoleFiller<TImage>::SetImage(typename TImage::Pointer image)
{
  this->Image = image;
  this->Output = TImage::New();
}

template <typename TImage>
void SmallHoleFiller<TImage>::SetHolePixel(typename TImage::PixelType pixel)
{
  this->HolePixel = pixel;
}

template <typename TImage>
typename TImage::Pointer SmallHoleFiller<TImage>::GetOutput()
{
  return this->Output;
}

template <typename TImage>
void SmallHoleFiller<TImage>::Fill()
{
  if(!this->Image)
    {
    std::cerr << "You must first call SetImage!" << std::endl;
    return;
    }

  // Initialize by setting the output image to the input image.
  DeepCopy<TImage>(this->Image, this->Output);
  unsigned int numberOfIterations = 0;
  while(HasEmptyPixels())
    {
    std::cout << "Iteration " << numberOfIterations << "..." << std::endl;
    Iterate();
    numberOfIterations++;
  
//    typedef  itk::ImageFileWriter< TImage  > WriterType;
//    typename WriterType::Pointer writer = WriterType::New();
//    std::stringstream ss;
//    ss << "/tmp/intermediate_" << numberOfIterations << ".mha";
//    writer->SetFileName(ss.str());
//    writer->SetInput(this->Output);
//    writer->Update();
    }
    
  std::cout << "Filling completed in " << numberOfIterations << " iterations." << std::endl;
}

template <typename TImage>
void SmallHoleFiller<TImage>::Iterate()
{
  // Make a copy of the current intermediate output. We use this to determine which pixels were holes at the beginning of this iteration.
  typename TImage::Pointer temporaryImage = TImage::New();
  DeepCopy<TImage>(this->Output, temporaryImage);
  
  // Traverse the image. When a pixel is encountered that is a hole, fill it with the average of its non-hole neighbors.
  // Do not mark pixels filled on this iteration as known. This will result in a less accurate filling, as it favors the colors
  // of pixels that occur earlier in the raster scan order.
  typename TImage::SizeType radius;
  radius.Fill(1);
  
  itk::NeighborhoodIterator<TImage> intermediateNeighborhoodIterator(radius, temporaryImage, temporaryImage->GetLargestPossibleRegion());
  
  // This iterator is used to set the output pixels.
  itk::ImageRegionIterator<TImage> outputIterator(this->Output, this->Output->GetLargestPossibleRegion());
    
  while(!intermediateNeighborhoodIterator.IsAtEnd())
    {
    if(intermediateNeighborhoodIterator.GetCenterPixel() == this->HolePixel) // We want to fill this pixel
      {
      typename TImage::PixelType pixelSum = itk::NumericTraits< typename TImage::PixelType >::Zero;
      // Loop over the 8-neighorhood
      unsigned int validPixels = 0;
      for(unsigned int i = 0; i < 27; i++)
	{
	if(i==13) // this is the center (current) pixel, so skip it
	  {
	  continue;
	  }
	bool isInBounds = false;
	typename TImage::PixelType currentPixel = intermediateNeighborhoodIterator.GetPixel(i, isInBounds);
	if(isInBounds)
	  {
	  if(currentPixel != this->HolePixel)
	    {
	    validPixels++;
	    pixelSum += currentPixel;
	    }
	  }
	} // end 8-connected neighbor for
	
      if(validPixels > 0)
	{
	//typename TImage::PixelType pixelAverage = static_cast<typename TImage::PixelType>(pixelSum / validPixels);
	typename TImage::PixelType pixelAverage = static_cast<typename TImage::PixelType>(pixelSum * (1.0/ validPixels)); // We multiply by the reciprocal because operator/ is not defined for all types.
	outputIterator.Set(pixelAverage);
	//std::cout << "Set " << outputIterator.GetIndex() << " to " << pixelAverage << std::endl;
	}
      } // end "fill this pixel?" conditional
      
    ++intermediateNeighborhoodIterator;
    ++outputIterator;
    } // end main traversal loop
    
}

template <typename TImage>
bool SmallHoleFiller<TImage>::HasEmptyPixels()
{
  itk::ImageRegionConstIterator<TImage> imageIterator(this->Output, this->Output->GetLargestPossibleRegion());
 
  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.Get() == this->HolePixel)
      {
      //std::cout << "Pixel " << imageIterator.GetIndex() << " is empty." << std::endl;
      return true;
      }
    ++imageIterator;
    }
  return false;
}

///////// This is not a member function /////////////
/** Copy the input to the output*/
template<typename TImage>
void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output)
{
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}
