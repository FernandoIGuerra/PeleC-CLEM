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
  amrex::ParmParse pp("prob");
  pp.query("pamb", PeleC::h_prob_parm_device->pamb);

  pp.query("T_cross", PeleC::h_prob_parm_device->T_cross);
  pp.query("u_cross", PeleC::h_prob_parm_device->u_cross);
  pp.query("Y_O2_cross", PeleC::h_prob_parm_device->Y_O2_cross);
  pp.query("Y_N2_cross", PeleC::h_prob_parm_device->Y_N2_cross);

  pp.query("T_jet", PeleC::h_prob_parm_device->T_jet);
  pp.query("v_jet", PeleC::h_prob_parm_device->v_jet);
  pp.query("jet_xc", PeleC::h_prob_parm_device->jet_xc);
  pp.query("jet_width", PeleC::h_prob_parm_device->jet_width);

  PeleC::h_prob_parm_device->jet_xlo =
    PeleC::h_prob_parm_device->jet_xc -
    0.5 * PeleC::h_prob_parm_device->jet_width;
  PeleC::h_prob_parm_device->jet_xhi =
    PeleC::h_prob_parm_device->jet_xc +
    0.5 * PeleC::h_prob_parm_device->jet_width;

  amrex::Print() << "CLEMTest2D: jet-in-crossflow reactor\n"
                 << "  crossflow: T = " << PeleC::h_prob_parm_device->T_cross
                 << " K, u = " << PeleC::h_prob_parm_device->u_cross
                 << " cm/s (air)\n"
                 << "  jet:       T = " << PeleC::h_prob_parm_device->T_jet
                 << " K, v = " << PeleC::h_prob_parm_device->v_jet
                 << " cm/s (pure H2), slot x in ["
                 << PeleC::h_prob_parm_device->jet_xlo << ", "
                 << PeleC::h_prob_parm_device->jet_xhi << "]\n";
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
