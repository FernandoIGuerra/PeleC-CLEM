#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("p_init", PeleC::h_prob_parm_device->p_init);
    pp.query("Y_H2_C", PeleC::h_prob_parm_device->Y_H2_C);
    pp.query("Y_O2_C", PeleC::h_prob_parm_device->Y_O2_C);
    pp.query("Y_N2_C", PeleC::h_prob_parm_device->Y_N2_C);
    pp.query("T_C", PeleC::h_prob_parm_device->T_C);

    pp.query("Y_H2_H", PeleC::h_prob_parm_device->Y_H2_H);
    pp.query("Y_O2_H", PeleC::h_prob_parm_device->Y_O2_H);
    pp.query("Y_N2_H", PeleC::h_prob_parm_device->Y_N2_H);
    pp.query("T_H", PeleC::h_prob_parm_device->T_H);

    pp.query("xliml", PeleC::h_prob_parm_device->xliml);
    pp.query("xlimh", PeleC::h_prob_parm_device->xlimh);
  }

  // Initial values
  PeleC::h_prob_parm_device->massfrac_C[H2_ID] = PeleC::h_prob_parm_device->Y_H2_C;
  PeleC::h_prob_parm_device->massfrac_C[O2_ID] = PeleC::h_prob_parm_device->Y_O2_C;
  PeleC::h_prob_parm_device->massfrac_C[N2_ID] = PeleC::h_prob_parm_device->Y_N2_C;

  PeleC::h_prob_parm_device->massfrac_H[H2_ID] = PeleC::h_prob_parm_device->Y_H2_H;
  PeleC::h_prob_parm_device->massfrac_H[O2_ID] = PeleC::h_prob_parm_device->Y_O2_H;
  PeleC::h_prob_parm_device->massfrac_H[N2_ID] = PeleC::h_prob_parm_device->Y_N2_H;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_init,
    PeleC::h_prob_parm_device->massfrac_C.begin(),
    PeleC::h_prob_parm_device->T_C, PeleC::h_prob_parm_device->rho_C,
    PeleC::h_prob_parm_device->e_C);

    eos.PYT2RE(
    PeleC::h_prob_parm_device->p_init,
    PeleC::h_prob_parm_device->massfrac_H.begin(),
    PeleC::h_prob_parm_device->T_H, PeleC::h_prob_parm_device->rho_H,
    PeleC::h_prob_parm_device->e_H);

}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
