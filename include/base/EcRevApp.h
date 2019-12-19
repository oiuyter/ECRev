#ifndef FLUXBCUDOT2_H
#define FLUXBCUDOT2_H

#include "IntegratedBC.h"

// Forward Declarations
class FluxBCudot2;
class Function;

template <>
InputParameters validParams<FluxBCudot2>();

/**
 * Implements Neumann BC where grad(u)=udot-something on the boundary.
 */

class FluxBCudot2 : public IntegratedBC

{
public:

  FluxBCudot2(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;


private:

  const VariableValue & _couple_var;
  const unsigned _couple_var_offjac;
  const VariableValue & _u_dot;
  const VariableValue & _du_dot_du;
  Real _K_L;
  Real _K_D;
  Real _G_S;
  Real _kf_1;
  Real _kb_1;
  Real _kf_2;
  Real _kb_2;
  Real _sigma;
  Real theta();
  Real gamma();
  Real diftheta();
  Real difgamma();
  Function & _func_1;
  Function & _func_2;
  Function & _func_3;
  Function & _func_4;
};

#endif
