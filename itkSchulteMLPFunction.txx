namespace itk
{

SchulteMLPFunction
::SchulteMLPFunction()
{
  // We operate a change of origin, u0 is always 0
  m_u0=0.;

  // Construct the constant part of R0 and R1 (equations 11 and 14)
  m_R0(0,0) = 1.;
  m_R0(1,0) = 0.;
  m_R0(1,1) = 1.;
  m_R1 = m_R0;

  // Transpose
  m_R0T = m_R0.GetTranspose();
  m_R1T = m_R1.GetTranspose();
}

void
SchulteMLPFunction
::Init(VectorType posIn, VectorType posOut, VectorType dirIn, VectorType dirOut)
{
  m_uOrigin = posIn[2];
  //m_IntForSigmaSqTheta0  = Functor::SchulteMLP::IntegralForSigmaSqTheta ::GetValue(m_u0);
  //m_IntForSigmaSqTTheta0 = Functor::SchulteMLP::IntegralForSigmaSqTTheta::GetValue(m_u0);
  //m_IntForSigmaSqT0      = Functor::SchulteMLP::IntegralForSigmaSqT     ::GetValue(m_u0);

  m_u2 = posOut[2]-m_uOrigin;
  m_IntForSigmaSqTheta2  = Functor::SchulteMLP::IntegralForSigmaSqTheta ::GetValue(m_u2);
  m_IntForSigmaSqTTheta2 = Functor::SchulteMLP::IntegralForSigmaSqTTheta::GetValue(m_u2);
  m_IntForSigmaSqT2      = Functor::SchulteMLP::IntegralForSigmaSqT     ::GetValue(m_u2);

  // Parameters vectors
  m_x0[0] = posIn[0];
  m_x0[1] = atan(dirIn[0]);  //dirIn[2] is implicitely 1.
  m_x2[0] = posOut[0];
  m_x2[1] = atan(dirOut[0]); //dirOut[2] is implicitely 1.

  m_y0[0] = posIn[1];
  m_y0[1] = atan(dirIn[1]);  //dirIn[2] is implicitely 1.
  m_y2[0] = posOut[1];
  m_y2[1] = atan(dirOut[1]); //dirOut[2] is implicitely 1.
}

void
SchulteMLPFunction
::Evaluate( const double u, double &x, double&y )
{
  const double u1 = u-m_uOrigin;

  // Finish constructing rotation matrices (equations 11 and 14)
  m_R0(0,1) = u1;
  m_R1(0,1) = m_u2-u1;
  m_R0T(1,0) = m_R0(0,1);
  m_R1T(1,0) = m_R1(0,1);

  // Constants used in both integrals
  const double intForSigmaSqTheta1  = Functor::SchulteMLP::IntegralForSigmaSqTheta ::GetValue(u1);
  const double intForSigmaSqTTheta1 = Functor::SchulteMLP::IntegralForSigmaSqTTheta::GetValue(u1);
  const double intForSigmaSqT1      = Functor::SchulteMLP::IntegralForSigmaSqT     ::GetValue(u1);

  // Construct Sigma1 (equations 6-9)
  itk::Matrix<double, 2, 2> Sigma1;
  Sigma1(1,1) = intForSigmaSqTheta1/* - m_IntForSigmaSqTheta0*/;
  Sigma1(0,1) = u1 * Sigma1(1,1) - intForSigmaSqTTheta1/* + m_IntForSigmaSqTTheta0*/;
  Sigma1(1,0) = Sigma1(0,1);
  Sigma1(0,0) = u1 * ( 2*Sigma1(0,1) - u1*Sigma1(1,1) ) + intForSigmaSqT1/* - m_IntForSigmaSqT0*/;
  Sigma1 *= Functor::SchulteMLP::ConstantPartOfIntegrals::GetValue(m_u0,u1);

  // Construct Sigma2 (equations 15-18)
  itk::Matrix<double, 2, 2> Sigma2;
  Sigma2(1,1) = m_IntForSigmaSqTheta2 - intForSigmaSqTheta1;
  Sigma2(0,1) = m_u2 * Sigma2(1,1) - m_IntForSigmaSqTTheta2 + intForSigmaSqTTheta1;
  Sigma2(1,0) = Sigma2(0,1);
  Sigma2(0,0) = m_u2 * ( 2*Sigma2(0,1) - m_u2*Sigma2(1,1) ) + m_IntForSigmaSqT2 - intForSigmaSqT1;
  Sigma2 *= Functor::SchulteMLP::ConstantPartOfIntegrals::GetValue(u1,m_u2);

  // x and y, equation 24
  // common calculations
  Sigma1 = Sigma1.GetInverse();
  Sigma2 = Sigma2.GetInverse();
  itk::Matrix<double, 2, 2> Sigma1Inv_R0 = Sigma1 * m_R0;
  itk::Matrix<double, 2, 2> R1T_Sigma2Inv = m_R1T * Sigma2;
  itk::Matrix<double, 2, 2> part(Sigma1 + R1T_Sigma2Inv * m_R1);
  part = part.GetInverse();

  // x
  itk::Vector<double, 2> xMLP;
  xMLP = part * (Sigma1Inv_R0 * m_x0 + R1T_Sigma2Inv * m_x2);
  x = xMLP[0];

  // y
  itk::Vector<double, 2> yMLP;
  yMLP = part * (Sigma1Inv_R0 * m_y0 + R1T_Sigma2Inv * m_y2);
  y = yMLP[0];
}

}
