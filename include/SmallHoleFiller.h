#ifndef SmallHoleFiller_H
#define SmallHoleFiller_H

template<typename TImage>
class SmallHoleFiller
{
public:
  // Constructor
  SmallHoleFiller();
  
  // Inputs
  void SetImage(typename TImage::Pointer image);
  void SetHolePixel(typename TImage::PixelType pixel);
  
  // Outputs
  typename TImage::Pointer GetOutput();
  
  // This is the main loop. It simply calls Iterate() until complete.
  void Fill();
  
  // This is the core functionality.
  void Iterate();
  
  // This function returns true if any of the Output pixels match the HolePixel. This indicates there is more work to be done.
  bool HasEmptyPixels();
  
private:
  // The input image.
  typename TImage::Pointer Image;
  
  // The intermediate and eventually output image.
  typename TImage::Pointer Output;
  
  // The value of a pixel to be considered a hole (to be filled).
  typename TImage::PixelType HolePixel;
  
  
};

// This function copies the data from 'input' to 'output'
template<typename TImage>
void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output);
  
  
#include "SmallHoleFiller.hxx"

#endif
