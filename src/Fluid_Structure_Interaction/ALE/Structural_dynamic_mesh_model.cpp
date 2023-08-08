#include <Structural_dynamic_mesh_model.h>

#include <Domaine_ALE.h>
#include <stdexcept>

#include <MGIS/Behaviour/Integrate.hxx>

Implemente_instanciable_sans_constructeur_ni_destructeur(Structural_dynamic_mesh_model,"Structural_dynamic_mesh_model",Interprete_geometrique_base);

Structural_dynamic_mesh_model::Structural_dynamic_mesh_model()
{

  l_ = "unset string" ;
  f_ = "unset string" ;
  h_ = "unset string" ;

  density_ = 0. ;
  inertialDamping_ = 0. ;
  youngModulus_ = 0. ;

  gridNStep = 0 ;
  gridTime = 0. ;
  gridDt = 0. ;
  isMassBuilt = false ;

}
Structural_dynamic_mesh_model::~Structural_dynamic_mesh_model()
{

}

Sortie& Structural_dynamic_mesh_model::printOn(Sortie& os) const
{
  return os;
}


Entree& Structural_dynamic_mesh_model::readOn(Entree& is)
{

}

Entree& Structural_dynamic_mesh_model::interpreter_(Entree& is)
{
  associer_domaine(is);
  Domaine_ALE& dom=ref_cast(Domaine_ALE, domaine());
  dom.reading_structural_dynamic_mesh_model(is);
  return is;
}

void Structural_dynamic_mesh_model::initMfrontBehaviour()
{

  if ( l_ == "unset string")
    {
      Cerr << "Error: Mfront library path undefined (Mfront_library keyword)" << finl ;
      exit() ;
    }

  if ( f_ == "unset string")
    {
      Cerr << "Error: Mfront model name undefined (Mfront_model_name keyword)" << finl ;
      exit() ;
    }

  if (h_ == "unset string")
    {
      switch (dimension)
        {
        case (2) :
          Cerr << "Use default Mfront hypothesis for 2D (PlaneStrain)" << finl ;
          h_ = "PlaneStrain" ;
          break ;
        case (3) :
          Cerr << "Use default Mfront hypothesis for 3D (Tridimensional)" << finl ;
          h_ = "Tridimensional" ;
          break ;
        }

    }
  hypothesis_= mgis::behaviour::fromString(h_) ;

  // dt/rdt
  mgisBehaviourData_.dt  = 0.;
  mgisBehaviourData_.rdt = &rdt_;

  mgisBehaviourData_.s0.material_properties = nullptr;
  mgisBehaviourData_.s1.material_properties = nullptr;

  load_behaviour_(StressMeasureKind::cauchy, TangentOperatorKind::none) ;
  if ( nbIvars_ > 0)
    {
      Cerr << "Wrong Mfront behaviour for structural dynamic mesh motion" << finl ;
      Cerr << "Only pure elastic materials are allowed" << finl ;
      exit();
    }
  if ( nbEvars_ > 0)
    {
      Cerr << "Wrong Mfront behaviour for structural dynamic mesh motion" << finl ;
      Cerr << "No external variable allowed" << finl ;
      exit();
    }

  // Check and fill material properties for Mfront behaviour
  matpSize_                                 = mgis::behaviour::getArraySize(mgisBehaviour_.mps, hypothesis_);
  double *s0MaterialProperties              = new double[matpSize_];
  double *s1MaterialProperties              = new double[matpSize_];
  mgisBehaviourData_.s0.material_properties = s0MaterialProperties;
  mgisBehaviourData_.s1.material_properties = s1MaterialProperties;
  for (auto mps : mgisBehaviour_.mps)
    {
      size_t vOffset = mgis::behaviour::getVariableOffset(mgisBehaviour_.mps, mps.name, hypothesis_);
      size_t vSize   = mgis::behaviour::getVariableSize(mps, hypothesis_);
      auto it        = matp_.find(mps.name);
      if (it != matp_.end())
        {
          if (vSize != (it->second).size())
            {
              Cerr << "Error with input for Mfront behaviour " << f_ << finl ;
              Cerr << "The number of scalar values must be " << vSize << " for property " << mps.name << finl ;
              exit();
            }
          for (size_t i = 0; i < vSize; ++i)
            {
              s0MaterialProperties[vOffset + i] = it->second[i];
              s1MaterialProperties[vOffset + i] = it->second[i];
            }
        }
      else
        {
          Cerr << "Error with input for Mfront behaviour " << f_ << finl ;
          Cerr << "Missing data for for property " << mps.name << finl ;
          exit();
        }
    }

}

void Structural_dynamic_mesh_model:: initDynamicMeshProblem(const int nsom, const int nelem, const MD_Vector& md)
{
  x.resize(nsom,dimension) ;
  u.resize(nsom,dimension) ;
  v.resize(nsom,dimension) ;
  vp.resize(nsom,dimension) ;
  a.resize(nsom,dimension) ;
  ff.resize(nsom,dimension) ;

  mass.resize(nsom) ;

  B0.resize(nelem,dimension*nbn_) ;
  Ft.resize(nelem, nSymSize_) ;
  Stress.resize(nelem, symSize_) ;

  massElem_.resize(nelem) ;

  // Only ff and mass need to be built as distributed arrays
  MD_Vector_tools::creer_tableau_distribue(md, ff) ;
  MD_Vector_tools::creer_tableau_distribue(md, mass) ;

  u=0 ;
  v=0 ;
  vp=0 ;
  a=0 ;
  ff=0 ;
  mass=0. ;

  Ft = 0. ;
  Stress = 0. ;
  massElem_ = 0. ;

  switch (dimension)
    {
    case (2) :
      nbn_ = 3 ;
      nSymSize_ = 4 ;
      symSize_ = 3 ;

      Eta_.resize(dimension,nbn_) ;
      Eta_(1,1) = -1. ;
      Eta_(1,2) =  1. ;
      Eta_(1,3) =  0. ;
      Eta_(2,1) = -1. ;
      Eta_(2,2) =  0. ;
      Eta_(2,3) =  1. ;

      break ;
    case (3) :
      nbn_ = 4 ;
      nSymSize_ = 9 ;
      symSize_ = 6 ;

      Eta_.resize(dimension,nbn_) ;
      Eta_(1,1) = -1. ;
      Eta_(1,2) =  1. ;
      Eta_(1,3) =  0. ;
      Eta_(1,4) =  0. ;
      Eta_(2,1) = -1. ;
      Eta_(2,2) =  0. ;
      Eta_(2,3) =  1. ;
      Eta_(2,4) =  0. ;
      Eta_(3,1) = -1. ;
      Eta_(3,2) =  0. ;
      Eta_(3,3) =  0. ;
      Eta_(3,4) =  1. ;

      break ;
    }

  xl_.resize(nbn_,dimension) ;
  ul_.resize(nbn_,dimension) ;
  fl_.resize(nbn_,dimension) ;

  B0l_.resize(dimension,nbn_) ;

}

void Structural_dynamic_mesh_model::computeInternalForces(double Vol, double Xlong, double E)
{

  // Integration with 1 single gauss point on triangular (dimension=2) or tetrahedral (dimension=3) finite element

  //  Compute jacobian and B matrices

  DoubleTab Jac(dimension, dimension) ;

  Jac=0 ;
  for (int k=0; k < nbn_; k++)
    {
      for (int j=0; j < dimension; j++)
        {
          for (int i=0; i < dimension ; i++)
            {
              Jac(j,i) += Eta_(j,k) * xl_(i,k) ;
            }
        }
    }

  double Det ;
  DoubleTab InvJac(dimension, dimension) ;
  double aux ;
  switch (dimension)
    {
    case (2) :
      Det = Jac(0,0) * Jac(1,1) - Jac(0,1) * Jac(1,0) ;
      aux = 1. / Det ;
      InvJac(0,0) =  Jac(1,1) * aux ;
      InvJac(0,1) = -Jac(0,1) * aux ;
      InvJac(1,0) = -Jac(1,0) * aux ;
      InvJac(1,1) =  Jac(0,0) * aux ;

      triangleSurfLength_(Vol, Xlong, Det) ;

      break ;
    case (3) :
      Det= Jac(0,0) * (Jac(1,1) * Jac(2,2) - Jac(1,2) * Jac(2,1))
           - Jac(0,1) * (Jac(1,0) * Jac(2,2) - Jac(1,2) * Jac(2,0))
           + Jac(0,2) * (Jac(1,0) * Jac(2,1) - Jac(1,1) * Jac(2,0)) ;
      aux = 1. / Det ;
      InvJac(0,0) =  (Jac(1,1) * Jac(2,2) - Jac(2,1) * Jac(1,2)) * aux ;
      InvJac(0,1) = -(Jac(0,1) * Jac(2,2) - Jac(2,1) * Jac(0,2)) * aux ;
      InvJac(0,2) =  (Jac(0,1) * Jac(1,2) - Jac(1,1) * Jac(0,2)) * aux ;
      InvJac(1,0) = -(Jac(1,0) * Jac(2,2) - Jac(2,0) * Jac(1,2)) * aux ;
      InvJac(1,1) =  (Jac(0,0) * Jac(2,2) - Jac(2,0) * Jac(0,2)) * aux ;
      InvJac(1,2) = -(Jac(0,0) * Jac(1,2) - Jac(1,0) * Jac(0,2)) * aux ;
      InvJac(2,0) =  (Jac(1,0) * Jac(2,1) - Jac(2,0) * Jac(1,1)) * aux ;
      InvJac(2,1) = -(Jac(0,0) * Jac(2,1) - Jac(2,0) * Jac(0,1)) * aux ;
      InvJac(2,2) =  (Jac(0,0) * Jac(1,1) - Jac(1,0) * Jac(0,1)) * aux ;

      tetrahedronVolLength_(Vol, Xlong, Det) ;

      break ;
    }

  DoubleTab B(dimension, nbn_) ;
  B=0 ;
  for (int k=0; k < dimension; k++)
    {
      for (int j=0; j < dimension; j++)
        {
          for (int i=0; i < nbn_ ; i++)
            {
              B(j,i) += InvJac(j,k) * Eta_(k,i) ;
            }
        }
    }

  if (gridNStep == 0) setB0_(B) ; // save B matrix on initial configuration

  //  Next transformation gradient

  double ftpdt[nSymSize_] ;
  for (int i=0; i<nSymSize_; i++) { ftpdt[i] = 0 ; }

#define _xikbjk(i,j) xl_(i,k)*B0l_(j,k)
  for (int k=0; k<nbn_; k++)
    {
      ftpdt[0] += _xikbjk(0,0); // F11=x1k*dNk/dx01
      ftpdt[1] += _xikbjk(1,1); // F22=x2k*dNk/dx02
      ftpdt[3] += _xikbjk(0,1); // F12=x1k*dNk/dx02
      ftpdt[4] += _xikbjk(1,0); // F21=x2k*dNk/dx01
      if (dimension == 3)
        {
          ftpdt[2] += _xikbjk(2,2); // F33=x3k*dNk/dx03
          ftpdt[5] += _xikbjk(0,2); // F13=x1k*dNk/dx03
          ftpdt[6] += _xikbjk(2,0); // F31=x3k*dNk/dx01
          ftpdt[7] += _xikbjk(1,2); // F23=x2k*dNk/dx03
          ftpdt[8] += _xikbjk(2,1); // F32=x3k*dNk/dx02
        }
    }

  //  Integrate material behavior

  double ft[nSymSize_] ;
  for (int i=0; i<nSymSize_; i++) { ft[i] = Ft(iel_,i) ; }

  double stress[symSize_] ;
  for (int i=0; i<symSize_; i++) { stress[i] = Stress(iel_,i) ; }

  double voidInternalVars[0] ;
  double voidExternalVars[0] ;
  double K[stiffnessMatrixMinSize_] ;
  K[0] = 0.; // Stiffness matrix not computed, see mfront/mgis convention

  integrate_behaviour_(stress, voidInternalVars, voidExternalVars, K, ft, ftpdt) ;

  // Update stress and transformation gradient
  for (int i=0; i<symSize_; i++) { Stress(iel_,i) = stress[i] ; }
  for (int i=0; i<nSymSize_; i++) { Ft(iel_,i) = ftpdt[i] ; }

  //  Compute local forces
  aux = 1./sqrt(2.) ;
  fl_ = 0 ;
  double cauchy[symSize_] ;
  switch (dimension)
    {
    case (2) :
      cauchy[0] = Stress[0] ;
      cauchy[1] = Stress[1] ;
      cauchy[2] = Stress[2] * aux ;

      for (int k=0; k<nbn_; k++)
        {
          fl_(0,k) += Vol * (B(0,k) * cauchy[0] + B(1,k) * cauchy[2]) ;
          fl_(1,k) += Vol * (B(0,k) * cauchy[2] + B(1,k) * cauchy[1]) ;
        }

      break ;
    case (3) :
      cauchy[0] = Stress[0] ;
      cauchy[1] = Stress[1] ;
      cauchy[2] = Stress[2] ;
      cauchy[3] = Stress[3] * aux ;
      cauchy[4] = Stress[4] * aux ;
      cauchy[5] = Stress[5] * aux ;

      for (int k=0; k<nbn_; k++)
        {
          fl_(0,k) += Vol * (B(0,k) * cauchy[0] + B(1,k) * cauchy[3] + B(2,k) * cauchy[5]) ;
          fl_(1,k) += Vol * (B(0,k) * cauchy[3] + B(1,k) * cauchy[1] + B(2,k) * cauchy[4]) ;
          fl_(2,k) += Vol * (B(0,k) * cauchy[5] + B(1,k) * cauchy[4] + B(2,k) * cauchy[2]) ;
        }

      break ;
    }

  // Temporary export of Young Modulus
  E = youngModulus_ ;

}

void Structural_dynamic_mesh_model::triangleSurfLength_(double aire, double xlong, const double Det)
{

  aire = Det / 2. ;

  DoubleVect v12(2) ;
  DoubleVect v13(2) ;
  DoubleVect v23(2) ;

  for (int i=0; i<dimension; i++)
    {
      v12[i] = xl_(i,1) - xl_(i,0) ;
      v13[i] = xl_(i,2) - xl_(i,0) ;
      v23[i] = xl_(i,2) - xl_(i,1) ;
    }

  DoubleVect w(3) ;

  w[0] = v12[0] * v12[0] + v12[1] * v12[1] ;
  w[1] = v13[0] * v13[0] + v13[1] * v13[1] ;
  w[2] = v23[0] * v23[0] + v23[1] * v23[1] ;

  double max_value = std::max({w[0], w[1], w[2]}) ;
  xlong = Det / sqrt(max_value) ;

}

void Structural_dynamic_mesh_model::tetrahedronVolLength_(double vol, double xlong, const double Det)
{

  vol = Det / 6. ;

  DoubleVect v21(3) ;
  DoubleVect v31(3) ;
  DoubleVect v41(3) ;
  DoubleVect v32(3) ;
  DoubleVect v42(3) ;

  for (int i=0; i<dimension; i++)
    {
      v21[i] = xl_(i,0) - xl_(i,1) ;
      v31[i] = xl_(i,0) - xl_(i,2) ;
      v41[i] = xl_(i,0) - xl_(i,3) ;
      v32[i] = xl_(i,1) - xl_(i,2) ;
      v42[i] = xl_(i,1) - xl_(i,3) ;
    }

  DoubleVect w(3), s(4) ;

  // face 1 2 3 : cross product 21^31
  w[0] = v21[1]*v31[2] - v21[2]*v31[1] ;
  w[1] = v21[2]*v31[0] - v21[0]*v31[2] ;
  w[2] = v21[0]*v31[1] - v21[1]*v31[0] ;
  s[0] = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] ;
  // face 1 2 4 : cross product 21^41
  w[0] = v21[1]*v41[2] - v21[2]*v41[1] ;
  w[1] = v21[2]*v41[0] - v21[0]*v41[2] ;
  w[2] = v21[0]*v41[1] - v21[1]*v41[0] ;
  s[1] = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] ;
  // face 1 3 4 : cross product 31^41
  w[0] = v31[1]*v41[2] - v31[2]*v41[1] ;
  w[1] = v31[2]*v41[0] - v31[0]*v41[2] ;
  w[2] = v31[0]*v41[1] - v31[1]*v41[0] ;
  s[2] = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] ;
  // face 2 3 4 : cross product 32^42
  w[0] = v32[1]*v42[2] - v32[2]*v42[1] ;
  w[1] = v32[2]*v42[0] - v32[0]*v42[2] ;
  w[2] = v32[0]*v42[1] - v32[1]*v42[0] ;
  s[3] = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] ;

  double max_value = std::max({s[0], s[1], s[2], s[3]}) ;
  xlong = Det / sqrt(max_value) ;

}

void Structural_dynamic_mesh_model::setB0_(const DoubleTab& B)
{

  // Save B operator at initial step to compute transformation gradient with respect to initial configuration

  int lb = nbn_ * dimension ;
  int pos = 0 ;

  for (int j=0; j < dimension; j++)
    {
      for (int i=0; i < nbn_ ; i++)
        {
          B0(iel_,pos) = B(j,i) ;
          pos++ ;
        }
    }

}

void Structural_dynamic_mesh_model::load_behaviour_(StressMeasureKind smk, TangentOperatorKind tok)
{
  if (loaded_)
    return;
  else
    loaded_ = true;

  auto getMgisStressMeasureKind = [](StressMeasureKind k)
  {
    switch (k)
      {
      case StressMeasureKind::cauchy:
        return mgis::behaviour::FiniteStrainBehaviourOptions::CAUCHY;
      case StressMeasureKind::pk1:
        return mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
      case StressMeasureKind::pk2:
        return mgis::behaviour::FiniteStrainBehaviourOptions::PK2;
      case StressMeasureKind::conjugate:
        return mgis::behaviour::FiniteStrainBehaviourOptions::CAUCHY;    // to avoid warning
      }
    return mgis::behaviour::FiniteStrainBehaviourOptions::CAUCHY;    // to avoid warning
  };

  auto getMgisTangentOperatorKind = [](TangentOperatorKind k)
  {
    switch (k)
      {
      case TangentOperatorKind::none:
        return mgis::behaviour::FiniteStrainBehaviourOptions::DSIG_DF;    // default mgis tangent operator kind
      case TangentOperatorKind::dCauchydF:
        return mgis::behaviour::FiniteStrainBehaviourOptions::DSIG_DF;
      case TangentOperatorKind::dPk2dE:
        return mgis::behaviour::FiniteStrainBehaviourOptions::DS_DEGL;
      case TangentOperatorKind::dPk1dF:
        return mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF;
      }
    return mgis::behaviour::FiniteStrainBehaviourOptions::DSIG_DF;    // to avoid warning
  };

  mgisBehaviour_ = mgis::behaviour::load(
  {getMgisStressMeasureKind(smk), getMgisTangentOperatorKind(tok)}, l_, f_, mgis::behaviour::fromString(h_));

  // internal variables
  this->nbIvars_ = mgis::behaviour::getArraySize(mgisBehaviour_.isvs, hypothesis_);

  // external variables
  this->nbEvars_ = mgis::behaviour::getArraySize(mgisBehaviour_.esvs, hypothesis_);

  // K
  this->KSize_ = mgis::behaviour::getArraySize(mgisBehaviour_.gradients, hypothesis_) *
                 mgis::behaviour::getArraySize(mgisBehaviour_.thermodynamic_forces, hypothesis_);
}

void Structural_dynamic_mesh_model::integrate_behaviour_(double *const stress,                     // stress tensor
                                                         double *const ivar,                       // state variables
                                                         const double *const evar,                 // external variables
                                                         double *const stiff,                      // stiffness matrix
                                                         const double *const gradientOrStrain0,    // deformation gradient or strain at the beginning of the time step
                                                         const double *const gradientOrStrain1     // deformation gradient or strain at the end of the time step
                                                        )
{
  mgisBehaviourData_.s0.gradients = const_cast<double *>(gradientOrStrain0);
  mgisBehaviourData_.s1.gradients = const_cast<double *>(gradientOrStrain1);

  mgisBehaviourData_.s0.thermodynamic_forces = stress;
  mgisBehaviourData_.s1.thermodynamic_forces = stress;

  mgisBehaviourData_.s0.internal_state_variables = ivar;
  mgisBehaviourData_.s1.internal_state_variables = ivar;

  mgisBehaviourData_.s0.external_state_variables = const_cast<double *>(evar);
  mgisBehaviourData_.s1.external_state_variables = const_cast<double *>(evar);

  mgisBehaviourData_.K = stiff;

  int ierr = mgis::behaviour::integrate(mgisBehaviourData_, mgisBehaviour_);

  if (ierr == -1)
    {
      // Message::error("Integration of the behavior failed");
      Cerr << "Error: Integration of the behaviour with Mfront failed" << finl ;
      exit() ;
    }
  else if (ierr == 0)
    {
      // Message::warning("Integration of the behavior unreliable");
      Cerr << "Warning: Integration of the behaviour with Mfront unreliable" << finl ;
    }
}

void Structural_dynamic_mesh_model::setLocalFields(const int elnodes[4], int elem)
{

  iel_ = elem ;

  nbn_ = 3 ;
  if (dimension == 3) nbn_ = 4 ;

  for (int j=0; j<dimension; j++)
    {
      for (int i=0; i<nbn_; i++)
        {
          int isom = elnodes[i] ;
          xl_(i,j) = x(j,isom) ;
          ul_(i,j) = u(j,isom) ;
        }
    }
}

void Structural_dynamic_mesh_model::setGlobalFields(const int elnodes[4])
{

  for (int j=0; j<dimension; j++)
    {
      for (int i=0; i<nbn_; i++)
        {
          int isom = elnodes[i] ;
          ff(j,isom) += fl_(i,j) ;
        }
    }
}

double Structural_dynamic_mesh_model::computeCriticalDt(double Vol, double Xlong, double E)
{

  // Current density
  double currDensity = massElem_(iel_) / Vol ;

  // Speed of sound
  double cSound = sqrt(E / currDensity) ;

  return Xlong / cSound ;

}
