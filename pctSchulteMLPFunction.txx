namespace pct
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
::Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut)
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
  m_x0[1] = vcl_atan(dirIn[0]);  //dirIn[2] is implicitely 1.
  m_x2[0] = posOut[0];
  m_x2[1] = vcl_atan(dirOut[0]); //dirOut[2] is implicitely 1.

  m_y0[0] = posIn[1];
  m_y0[1] = vcl_atan(dirIn[1]);  //dirIn[2] is implicitely 1.
  m_y2[0] = posOut[1];
  m_y2[1] = vcl_atan(dirOut[1]); //dirOut[2] is implicitely 1.
}

void
SchulteMLPFunction
::Evaluate( const double u, double &x, double&y )
{
#ifdef MLP_TIMING
  m_EvaluateProbe1.Start();
#endif
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
  m_Sigma1(1,1) = intForSigmaSqTheta1/* - m_IntForSigmaSqTheta0*/;
  m_Sigma1(0,1) = u1 * m_Sigma1(1,1) - intForSigmaSqTTheta1/* + m_IntForSigmaSqTTheta0*/;
  m_Sigma1(1,0) = m_Sigma1(0,1);
  m_Sigma1(0,0) = u1 * ( 2*m_Sigma1(0,1) - u1*m_Sigma1(1,1) ) + intForSigmaSqT1/* - m_IntForSigmaSqT0*/;
  m_Sigma1 *= Functor::SchulteMLP::ConstantPartOfIntegrals::GetValue(m_u0,u1);

  // Construct Sigma2 (equations 15-18)
  m_Sigma2(1,1) = m_IntForSigmaSqTheta2 - intForSigmaSqTheta1;
  m_Sigma2(0,1) = m_u2 * m_Sigma2(1,1) - m_IntForSigmaSqTTheta2 + intForSigmaSqTTheta1;
  m_Sigma2(1,0) = m_Sigma2(0,1);
  m_Sigma2(0,0) = m_u2 * ( 2*m_Sigma2(0,1) - m_u2*m_Sigma2(1,1) ) + m_IntForSigmaSqT2 - intForSigmaSqT1;
  m_Sigma2 *= Functor::SchulteMLP::ConstantPartOfIntegrals::GetValue(u1,m_u2);

#ifdef MLP_TIMING
  m_EvaluateProbe1.Stop();
  m_EvaluateProbe2.Start();
#endif

  // x and y, equation 24
  // common calculations
  InverseMatrix(m_Sigma1);
  InverseMatrix(m_Sigma2);
  itk::Matrix<double, 2, 2> Sigma1Inv_R0 = m_Sigma1 * m_R0;
  itk::Matrix<double, 2, 2> R1T_Sigma2Inv = m_R1T * m_Sigma2;
  itk::Matrix<double, 2, 2> part(m_Sigma1 + R1T_Sigma2Inv * m_R1);
  InverseMatrix(part);

  // x
  itk::Vector<double, 2> xMLP;
  xMLP = part * (Sigma1Inv_R0 * m_x0 + R1T_Sigma2Inv * m_x2);
  x = xMLP[0];

  // y
  itk::Vector<double, 2> yMLP;
  yMLP = part * (Sigma1Inv_R0 * m_y0 + R1T_Sigma2Inv * m_y2);
  y = yMLP[0];

#ifdef MLP_TIMING
  m_EvaluateProbe2.Stop();
#endif
}

void
SchulteMLPFunction
::EvaluateError( const double u, itk::Matrix<double, 2, 2> &error )
{
  double x, y;
  Evaluate(u,x,y);
  error = m_Sigma1 + m_R1T * m_Sigma2 * m_R1;
  InverseMatrix(error);
  error *= 2.;
}

#ifdef MLP_TIMING
void
SchulteMLPFunction
::PrintTiming(std::ostream& os)
{
  os << "SchulteMLPFunction timing:" << std::endl;
  os << "  EvaluateProbe1: " << m_EvaluateProbe1.GetTotal()
     << ' ' << m_EvaluateProbe1.GetUnit() << std::endl;
  os << "  EvaluateProbe2: " << m_EvaluateProbe2.GetTotal()
     << ' ' << m_EvaluateProbe2.GetUnit() << std::endl;
}
#endif

void
SchulteMLPFunction
::InverseMatrix(itk::Matrix<double, 2, 2> &mat)
{
  double det = 1. / ( mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0) );
  std::swap( mat(0,0), mat(1,1) );
  mat(1,0) *= -1.;
  mat(0,1) *= -1.;
  mat *= det;
}

}
