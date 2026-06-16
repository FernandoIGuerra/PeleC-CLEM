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

    pp.query("Y_H2_low", PeleC::h_prob_parm_device->Y_H2_low);
    pp.query("Y_O2_low", PeleC::h_prob_parm_device->Y_O2_low);
    pp.query("Y_N2_low", PeleC::h_prob_parm_device->Y_N2_low);
    pp.query("T_low", PeleC::h_prob_parm_device->T_low);

    pp.query("Y_H2_high", PeleC::h_prob_parm_device->Y_H2_high);
    pp.query("Y_O2_high", PeleC::h_prob_parm_device->Y_O2_high);
    pp.query("Y_N2_high", PeleC::h_prob_parm_device->Y_N2_high);
    pp.query("T_high", PeleC::h_prob_parm_device->T_high);

    pp.query("yliml", PeleC::h_prob_parm_device->yliml);
    pp.query("ylimh", PeleC::h_prob_parm_device->ylimh);

    pp.query("velx_low", PeleC::h_prob_parm_device->velx_low);
    pp.query("vely_low", PeleC::h_prob_parm_device->vely_low);
    pp.query("velx_high", PeleC::h_prob_parm_device->velx_high);
    pp.query("vely_high", PeleC::h_prob_parm_device->vely_high);
  }

  // Initial values
  PeleC::h_prob_parm_device->massfrac_low[H2_ID] = PeleC::h_prob_parm_device->Y_H2_low;
  PeleC::h_prob_parm_device->massfrac_low[O2_ID] = PeleC::h_prob_parm_device->Y_O2_low;
  PeleC::h_prob_parm_device->massfrac_low[N2_ID] = PeleC::h_prob_parm_device->Y_N2_low;

  PeleC::h_prob_parm_device->massfrac_high[H2_ID] = PeleC::h_prob_parm_device->Y_H2_high;
  PeleC::h_prob_parm_device->massfrac_high[O2_ID] = PeleC::h_prob_parm_device->Y_O2_high;
  PeleC::h_prob_parm_device->massfrac_high[N2_ID] = PeleC::h_prob_parm_device->Y_N2_high;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_init,
    PeleC::h_prob_parm_device->massfrac_low.begin(),
    PeleC::h_prob_parm_device->T_low, PeleC::h_prob_parm_device->rho_low,
    PeleC::h_prob_parm_device->e_low);

    eos.PYT2RE(
    PeleC::h_prob_parm_device->p_init,
    PeleC::h_prob_parm_device->massfrac_high.begin(),
    PeleC::h_prob_parm_device->T_high, PeleC::h_prob_parm_device->rho_high,
    PeleC::h_prob_parm_device->e_high);

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
