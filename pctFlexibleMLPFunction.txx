namespace pct
{

FlexibleMLPFunction
::FlexibleMLPFunction()
{
  // We operate a change of origin, u0 is always 0
  m_u0=0.;
  m_ScalarTest = -1.;
  m_CanBeVectorised = true;
}

void
FlexibleMLPFunction
::InitUncertain(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut, double dEntry, double dExit, double m_TrackerResolution, double m_TrackerPairSpacing, double m_MaterialBudget)
{
  itkGenericExceptionMacro("Method InitUncertain not implemented for this derived class FlexibleMLPFunction.");
}

void
FlexibleMLPFunction
::Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut)
{
  itkGenericExceptionMacro("This version of the Init method not implemented for derived class PolynomialMLPFunction.");
}


void
FlexibleMLPFunction
::Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut, double eIn, double eOut)
{
  m_uOrigin = posIn[2];
  m_u2 = posOut[2]-m_uOrigin;

  // Parameters vectors
  m_x0[0] = posIn[0];
  m_x0[1] = std::atan(dirIn[0]);  //dirIn[2] is implicitely 1.
  m_x2[0] = posOut[0];
  m_x2[1] = std::atan(dirOut[0]); //dirOut[2] is implicitely 1.

  m_y0[0] = posIn[1];
  m_y0[1] = std::atan(dirIn[1]);  //dirIn[2] is implicitely 1.
  m_y2[0] = posOut[1];
  m_y2[1] = std::atan(dirOut[1]); //dirOut[2] is implicitely 1.

  Functor::FlexibleMLP::GetLinearCoefficients( m_ab, eIn, eOut, m_u2 );

  const double A = Functor::FlexibleMLP::FactorsABCD::GetA(m_u2, m_ab);
  const double B = Functor::FlexibleMLP::FactorsABCD::GetB(m_u2, m_ab);
  const double C = Functor::FlexibleMLP::FactorsABCD::GetC(m_u2, m_ab);
  const double D = Functor::FlexibleMLP::FactorsABCD::GetD(m_u2, m_ab);

  Functor::FlexibleMLP::CoefficientsC::GetValue(m_c_x, m_u2, m_x0, m_x2, A, B, C, D);
  Functor::FlexibleMLP::CoefficientsC::GetValue(m_c_y, m_u2, m_y0, m_y2, A, B, C, D);

  double bOvera = m_ab[1]/m_ab[0];
  double bOvera2 = m_ab[1]/m_ab[0]/m_ab[0];
  double b2Overa3 = bOvera2*m_ab[1]/m_ab[0];
  m_dm_x[0] = m_x0[0];
  m_dm_x[3] = m_c_x[0]*bOvera2 - m_c_x[1]*b2Overa3;
  m_dm_x[1] = m_x0[1]*bOvera - m_dm_x[3];
  m_dm_x[2] = 0.5*m_c_x[1]*b2Overa3;

  m_dm_y[0] = m_y0[0];
  m_dm_y[3] = m_c_y[0]*bOvera2 - m_c_y[1]*b2Overa3;
  m_dm_y[1] = m_y0[1]*bOvera - m_dm_y[3];
  m_dm_y[2] = 0.5*m_c_y[1]*b2Overa3;

}

void
FlexibleMLPFunction
::Evaluate( const double u, double &x, double&y, double &dx, double&dy )
{
  itkGenericExceptionMacro("Method Evaluate not implemented for this derived class FlexibleMLPFunction.");
}

// vectorised version
void
FlexibleMLPFunction
::Evaluate( std::vector<double> u, std::vector<double> &x, std::vector<double> &y )
{
  // shift so u starts at 0 and scale by a/b to get u_tilde
  for(auto& element : u)
  {
    element -= m_uOrigin;
    element *= (m_ab[0]/m_ab[1]);
  }

  // uLog = log(u_tilde + 1)
  std::vector<double> uLog;
  uLog.reserve(u.size());
  std::copy(u.begin(),u.end(),std::back_inserter(uLog));
  for(auto& element : uLog)
  {
    element += 1;
    element = std::log(element);
  }

  #ifdef MLP_TIMING
    m_EvaluateProbe1.Start();
  #endif

  /* terms order 0, 1, 2 */
  std::fill(x.begin(), x.end(), m_dm_x[2]);
  std::transform(x.begin(), x.end(), u.begin(), x.begin(), std::multiplies<double>() );
  std::transform(x.begin(), x.end(), x.begin(), bind2nd(std::plus<double>(), m_dm_x[1]));
  std::transform(x.begin(), x.end(), u.begin(), x.begin(), std::multiplies<double>() );
  std::transform(x.begin(), x.end(), x.begin(), bind2nd(std::plus<double>(), m_dm_x[0]));

  std::fill(y.begin(), y.end(), m_dm_y[2]);
  std::transform(y.begin(), y.end(), u.begin(), y.begin(), std::multiplies<double>() );
  std::transform(y.begin(), y.end(), y.begin(), bind2nd(std::plus<double>(), m_dm_y[1]));
  std::transform(y.begin(), y.end(), u.begin(), y.begin(), std::multiplies<double>() );
  std::transform(y.begin(), y.end(), y.begin(), bind2nd(std::plus<double>(), m_dm_y[0]));

  /* logarithmic term */
  std::transform(u.begin(), u.end(), u.begin(), bind2nd(std::plus<double>(), 1) ); // u => u+1
  std::transform(u.begin(), u.end(), uLog.begin(), uLog.begin(), std::multiplies<double>() ); // uLog = (u+1)*log(u+1)
  std::transform(uLog.begin(), uLog.end(), uLog.begin(), bind2nd(std::multiplies<double>(), m_dm_x[3]));
  std::transform(x.begin(), x.end(), uLog.begin(), x.begin(), std::plus<double>() );
  std::transform(uLog.begin(), uLog.end(), uLog.begin(), bind2nd(std::multiplies<double>(), m_dm_y[3]/m_dm_x[3]));
  std::transform(y.begin(), y.end(), uLog.begin(), y.begin(), std::plus<double>() );


#ifdef MLP_TIMING
  m_EvaluateProbe1.Stop();
#endif
}

void
FlexibleMLPFunction
::EvaluateError( const double u, itk::Matrix<double, 2, 2> &error )
{
  itkGenericExceptionMacro("The method PolynomialMLPFunction::EvaluateError is not implemented at the moment");
}

#ifdef MLP_TIMING
void
FlexibleMLPFunction
::PrintTiming(std::ostream& os)
{
  os << "PolynomialMLPFunction timing:" << std::endl;
  os << "  EvaluateProbe1: " << m_EvaluateProbe1.GetTotal()
     << ' ' << m_EvaluateProbe1.GetUnit() << std::endl;
  // os << "  EvaluateProbe2: " << m_EvaluateProbe2.GetTotal()
  //    << ' ' << m_EvaluateProbe2.GetUnit() << std::endl;
}
#endif


}
