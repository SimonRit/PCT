namespace pct
{

template < class TCoordRep >
void
ThirdOrderPolynomialMLPFunction<TCoordRep>
::Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut)
{
  // Parameters of the 3rd order polynomial. The function goes from
  // z=0 to z=lastSliceZ and an offset in the z direction is stored.
  zoffset = posIn[2];
  ax = posIn[0];
  ay = posIn[1];
  const TCoordRep x  = posOut[0];
  const TCoordRep y  = posOut[1];
  bx = dirIn[0];
  by = dirIn[1];
  const TCoordRep xd = dirOut[0];
  const TCoordRep yd = dirOut[1];
  const TCoordRep lastSliceZ = (posOut[2]-posIn[2]);
  const TCoordRep invzsq = 1./(lastSliceZ*lastSliceZ);
  cx = invzsq * ( 3*x - lastSliceZ*xd - 3*ax - 2*bx*lastSliceZ );
  cy = invzsq * ( 3*y - lastSliceZ*yd - 3*ay - 2*by*lastSliceZ );
  const TCoordRep inv3zsq = invzsq/3.;
  dx = inv3zsq * ( xd - bx - 2*cx*lastSliceZ );
  dy = inv3zsq * ( yd - by - 2*cy*lastSliceZ );
}

template < class TCoordRep >
void
ThirdOrderPolynomialMLPFunction<TCoordRep>
::Evaluate( const TCoordRep z, TCoordRep &x, TCoordRep&y )
{
  const TCoordRep zz = (z-zoffset);
  x = ax+zz*(bx+zz*(cx+zz*dx));
  y = ay+zz*(by+zz*(cy+zz*dy));
}

template < class TCoordRep >
void
ThirdOrderPolynomialMLPFunction<TCoordRep>
::Evaluate( std::vector<double> u, std::vector<double> &x, std::vector<double> &y )
{
  std::cout << "Vectorised version of Evaluate method not implemented for derived class SchulteMLPFunction." << std::endl;
}

}
