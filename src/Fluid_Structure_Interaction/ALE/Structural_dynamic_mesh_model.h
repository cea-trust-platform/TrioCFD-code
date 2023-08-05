/****************************************************************************/
#ifndef Structural_dynamic_mesh_model_included
#define Structural_dynamic_mesh_model_included

#include <Interprete_geometrique_base.h>
#include <stdexcept>

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourDataView.hxx>

// this enum defines the stress measure that the behavior is expected to return when calling 'integrate'
//   cauchy     -> cauchy stress tensor
//   pk1        -> first Piola-Kirchhoff tensor (aka Boussinesq tensor).
//                 This one is non-symmetric. It is such that P=pk1:dFdt, and Div(pk1^t)=f
//   pk2        -> second Piola-Kirchhoff tensor
//   conjugate  -> stress tensor conjugate with the input strain tensor, i.e. such that P=s:dedt
//                 where 'P' is the mechanical power, 's' is the stress tensor and 'dedt' the rate of the input strain tensor
enum class StressMeasureKind { cauchy, pk1, pk2, conjugate };

// this enum defines the tangent operator that the behavior is expected to return when calling 'integrate'
//   none      -> the tangent operator is not computed
//   dCauchydF -> derivative of the Cauchy stress tensor with respect to the gradient of the material deformation gradient
//   dPk2dE    -> derivative of the 2nd Piola-Kirchhoff stress tensor with respect to the Green-Lagrange strain tensor
//   dPk1dF    -> derivative of the 1st Piola-Kirchhoff stress tensor with respect to the gradient of the material deformation gradient
enum class TangentOperatorKind { none, dCauchydF, dPk2dE, dPk1dF };

class Structural_dynamic_mesh_model : public Interprete_geometrique_base
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Structural_dynamic_mesh_model);

public:
  Structural_dynamic_mesh_model();
  ~Structural_dynamic_mesh_model();
  Entree& interpreter_(Entree&) override;

  // Mfront behaviour
  // ----------------

  inline void setMfrontLibraryPath(const std::string& l) ;
  inline void setMfrontModelName(const std::string& h) ;

  // End Mfront behaviour
  // --------------------

  inline void setElasticityModulus(const double& E) ;
  inline void setPoissonRatio(const double& nu) ;
  inline void setInertialDamping(const double& D) ;
  inline void setDtSafetyCoefficient(const double& c) ;
  inline int getMinNumberOfMaterialProperties() const ;

  void initMfrontBehaviour() ;
  void initDynamicMeshProblem(const int nsom, const int nelem) ;

  void computeInternalForces() ;

protected:

  // Mfront behaviour
  // ----------------

  mgis::behaviour::Behaviour mgisBehaviour_;
  mgis::behaviour::BehaviourDataView mgisBehaviourData_;

  static constexpr int stiffnessMatrixMinSize_ = 1 + mgis::behaviour::Behaviour::nopts;

  void load_behaviour_(StressMeasureKind smk, TangentOperatorKind tok);


  std::string l_;    // behavior library path
  std::string f_;    // name of the behavior in the library
  std::string h_;    // modelling hypothesis (tridimentional, plane strain, plane stress, ...)
  /*
  Valid values for h_ are:
   - `AxisymmetricalGeneralisedPlaneStrain`
   - `AxisymmetricalGeneralisedPlaneStress`
   - `Axisymmetrical`
   - `PlaneStress`
   - `PlaneStrain`
   - `GeneralisedPlaneStrain`
   - `Tridimensional`
  */

  int minNumberOfMaterialProperties_ = 2 ;

  double elasticityModulus_ = 0;
  double poissonRatio_ = 0;

  double rdt_ = 1.;    // this is used by mgis

  // numerotation, mfront convention:
  // for symetric tensors:
  //  - 2D: [t00, t11, t22, sqrt(2)*t01]
  //  - 3D: [t00, t11, t22, sqrt(2)*t01, sqrt(2)*t02, sqrt(2)*t12]
  // for non-symetric tensors:
  //  - 2D: [t00, t11, t22, t01, t10]
  //  - 3D: [t00, t11, t22, t01, t10, t02, t20, t12, t21]
  void integrate_behaviour_(double *const stress,                     // stress tensor
                            double *const ivar,                       // state variables
                            const double *const evar,                 // external variables
                            double *const stiff,                      // stiffness matrix
                            const double *const gradientOrStrain0,    // deformation gradient or strain at the beginning of the time step
                            const double *const gradientOrStrain1     // deformation gradient or strain at the end of the time step
                           );

  // End Mfront behaviour
  // --------------------

  // Parameters and global variables
  // -------------------------------

  double inertialDamping_ = 0;
  double dtSafetyCoefficient_ = 0.8 ;

  int nstep_ ; // number of steps performed in structural dynamics integration
  double time_ ; // internal time counter

  int nbn_ ; // number of vertex nodes
  int nSymSize_ ; // length of non-symetric tensor in Voigt notation
  int symSize_ ; // length of symetric tensor in Voigt notation
  DoubleTab Eta_ ; // Shape function derivatives

  // Local variables for internal forces computation
  // -----------------------------------------------

  int iel_ ; // id of current element
  DoubleTab xl_ ; // local current coordinates
  DoubleTab ul_ ; // local nodal displacements
  DoubleTab fl_ ; // local nodal forces

  DoubleTab B0l_ ; // local B matrix on initial configuration

  // Global vectors and arrays for dynamic time integration
  // ------------------------------------------------------

  DoubleTab B0_ ;
  DoubleTab Ft_ ;
  DoubleTab Stress_ ;

  DoubleTab x_ ;
  DoubleTab u_ ;
  DoubleTab v_ ;
  DoubleTab vp_ ;
  DoubleTab a_ ;
  DoubleTab ff_ ;

  // Functions
  // ---------

  void triangleSurfLength_(double aire, double xlong, const double Det);
  void tetrahedronVolLength_(double vol, double xlong, const double Det);
  void setB0_(const DoubleTab& B) ;

private:

  // Mfront behaviour
  // ----------------

  int nbIvars_;
  int nbEvars_;
  int matpSize_;
  int KSize_;
  mgis::behaviour::Hypothesis hypothesis_;

  bool loaded_ = false;

  // End Mfront behaviour
  // --------------------

};

inline void Structural_dynamic_mesh_model::setMfrontLibraryPath(const std::string& l)
{
  l_ = l ;
}

inline void Structural_dynamic_mesh_model::setMfrontModelName(const std::string& h)
{
  h_ = h ;
}

inline void Structural_dynamic_mesh_model::setElasticityModulus(const double& E)
{
  elasticityModulus_ = E ;
}

inline void Structural_dynamic_mesh_model::setPoissonRatio(const double& nu)
{
  poissonRatio_ = nu ;
}

inline void Structural_dynamic_mesh_model::setInertialDamping(const double& D)
{
  inertialDamping_ = D ;
}

inline void Structural_dynamic_mesh_model::setDtSafetyCoefficient(const double& c)
{
  dtSafetyCoefficient_ = c ;
}

inline int Structural_dynamic_mesh_model::getMinNumberOfMaterialProperties() const
{
  return minNumberOfMaterialProperties_ ;
}

#endif
