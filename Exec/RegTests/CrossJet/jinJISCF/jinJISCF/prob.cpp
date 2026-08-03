#include "prob.H"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

namespace {
struct JetBoundaryProbeState
{
  bool enabled = false;
  bool append = false;
  bool header_written = false;
  int write_int = 100;
  std::string file = "jet_boundary_probe.csv";
  amrex::Real target_pressure = std::numeric_limits<amrex::Real>::quiet_NaN();
  amrex::Real target_temperature = 300.0;
  amrex::Real target_mass_flow_y = 0.109;
  amrex::Real target_mach = 1.02;
  amrex::Real target_velocity = 0.0;
};

JetBoundaryProbeState g_jet_probe;

void
write_jet_probe_row(
  const int step,
  const amrex::Real time,
  const amrex::Real sample_count,
  const amrex::Real jet_area,
  const amrex::Real avg_pressure,
  const amrex::Real avg_temperature,
  const amrex::Real avg_density,
  const amrex::Real avg_u,
  const amrex::Real avg_v,
  const amrex::Real avg_w,
  const amrex::Real avg_magvel,
  const amrex::Real mass_flow_y,
  const amrex::Real h2_mass_flow_y,
  const amrex::Real momentum_flux_y,
  const amrex::Real avg_momentum_flux_y,
  const amrex::Real crossflow_momentum_flux,
  const amrex::Real effective_momentum_ratio,
  const ProbParmDevice& prob_parm)
{
  std::ios_base::openmode mode = std::ios::out;
  if (g_jet_probe.header_written || g_jet_probe.append) {
    mode |= std::ios::app;
  }

  std::ofstream ofs(g_jet_probe.file, mode);
  ofs << std::setprecision(16);

  const bool need_header =
    !g_jet_probe.header_written &&
    (!g_jet_probe.append || ofs.tellp() == std::streampos(0));
  if (need_header) {
    ofs << "step,time,n_cells,centx,centz,d_jet,avg_pressure,"
           "avg_temperature,avg_density,avg_x_velocity,avg_y_velocity,"
           "avg_z_velocity,avg_mag_velocity,jet_area,mass_flow_y,"
           "h2_mass_flow_y,momentum_flux_y,avg_momentum_flux_y,"
           "crossflow_momentum_flux,effective_momentum_ratio,"
           "target_global_equivalence_ratio,"
           "target_jet_to_freestream_momentum_ratio,"
           "target_pressure,target_temperature,target_mass_flow_y,"
           "target_mach_number,"
           "target_x_velocity,target_y_velocity,target_z_velocity\n";
  }

  ofs << step << ',' << time << ',' << sample_count << ','
      << prob_parm.centx << ',' << prob_parm.centz << ','
      << 2.0 * prob_parm.r_hole << ',' << avg_pressure << ','
      << avg_temperature << ',' << avg_density << ',' << avg_u << ','
      << avg_v << ',' << avg_w << ',' << avg_magvel << ',' << jet_area
      << ',' << mass_flow_y << ',' << h2_mass_flow_y << ','
      << momentum_flux_y << ',' << avg_momentum_flux_y << ','
      << crossflow_momentum_flux << ',' << effective_momentum_ratio
      << ',' << prob_parm.global_equivalence_ratio << ','
      << prob_parm.jet_to_freestream_momentum_ratio << ','
      << g_jet_probe.target_pressure << ','
      << g_jet_probe.target_temperature << ','
      << g_jet_probe.target_mass_flow_y << ','
      << g_jet_probe.target_mach << ',' << 0.0 << ','
      << g_jet_probe.target_velocity << ',' << 0.0 << '\n';

  g_jet_probe.header_written = true;
}
} // namespace

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
  pp.query("Pres_amb", PeleC::h_prob_parm_device->Pres_amb);
  pp.query("Temp_amb", PeleC::h_prob_parm_device->Temp_amb);
  pp.query("Yox_amb", PeleC::h_prob_parm_device->Yox_amb);
  pp.query("Yox_channel", PeleC::h_prob_parm_device->Yox_channel);
  pp.query("Mach_channel", PeleC::h_prob_parm_device->Mach_channel);
  pp.query("Pres_jet", PeleC::h_prob_parm_device->Pres_jet);
  pp.query("temp_jet", PeleC::h_prob_parm_device->temp_jet);
  pp.query("Yfuel_jet", PeleC::h_prob_parm_device->Yfuel_jet);
  pp.query("Mach_jet", PeleC::h_prob_parm_device->Mach_jet);
  pp.query("vel_jet", PeleC::h_prob_parm_device->vel_jet);
  pp.query("massflow_jet", PeleC::h_prob_parm_device->massflow_jet);
  pp.query(
    "global_equivalence_ratio",
    PeleC::h_prob_parm_device->global_equivalence_ratio);
  pp.query(
    "jet_to_freestream_momentum_ratio",
    PeleC::h_prob_parm_device->jet_to_freestream_momentum_ratio);
  pp.query("inject_fuel", PeleC::h_prob_parm_device->inject_fuel);
  pp.query("centx", PeleC::h_prob_parm_device->centx);
  pp.query("centz", PeleC::h_prob_parm_device->centz);
  pp.query("r_hole", PeleC::h_prob_parm_device->r_hole);

  g_jet_probe.target_pressure = PeleC::h_prob_parm_device->Pres_jet;
  g_jet_probe.target_temperature = PeleC::h_prob_parm_device->temp_jet;
  g_jet_probe.target_mass_flow_y =
    PeleC::h_prob_parm_device->massflow_jet;
  g_jet_probe.target_mach = PeleC::h_prob_parm_device->Mach_jet;
  g_jet_probe.target_velocity = PeleC::h_prob_parm_device->vel_jet;
  // Reference velocity only. The actual BC velocity is computed in prob.H
  // from massflow_jet and the density implied by the interior pressure.
  pp.query("jet_probe_enable", g_jet_probe.enabled);
  pp.query("jet_probe_append", g_jet_probe.append);
  pp.query("jet_probe_int", g_jet_probe.write_int);
  pp.query("jet_probe_file", g_jet_probe.file);
  pp.query("jet_probe_target_pressure", g_jet_probe.target_pressure);
  pp.query("jet_probe_target_temperature", g_jet_probe.target_temperature);
  pp.query("jet_probe_target_mass_flow_y", g_jet_probe.target_mass_flow_y);
  pp.query("jet_probe_target_mach", g_jet_probe.target_mach);
  pp.query("jet_probe_target_velocity", g_jet_probe.target_velocity);
}

void
PeleC::problem_post_timestep()
{
  if (
    level != 0 || !g_jet_probe.enabled || g_jet_probe.write_int <= 0 ||
    PeleC::h_prob_parm_device->inject_fuel == 0) {
    return;
  }

  const int step = parent->levelSteps(0);
  if (step % g_jet_probe.write_int != 0) {
    return;
  }

  const amrex::Real time = state[State_Type].curTime();
  const amrex::MultiFab& S_new = get_new_data(State_Type);
  auto pressure_mf = derive("pressure", time, 0);

  constexpr int nvals = 12;
  amrex::Gpu::DeviceVector<amrex::Real> d_sums(nvals, 0.0);
  auto* sums = d_sums.data();

  const amrex::Geometry& geom0 = geom;
  amrex::Box ylo_box(geom0.Domain());
  ylo_box.setRange(1, geom0.Domain().smallEnd(1), 1);
  const auto prob_lo = geom0.ProbLoArray();
  const auto dx = geom0.CellSizeArray();
  const int ilo = geom0.Domain().smallEnd(0);
  const int klo = geom0.Domain().smallEnd(2);

  const amrex::Real hole_cx = PeleC::h_prob_parm_device->centx;
  const amrex::Real hole_cz = PeleC::h_prob_parm_device->centz;
  const amrex::Real hole_radius = PeleC::h_prob_parm_device->r_hole;
  const amrex::Real hole_radius2 = hole_radius * hole_radius;
  const int inject_fuel = PeleC::h_prob_parm_device->inject_fuel;

  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box bx = mfi.validbox() & ylo_box;
    if (!bx.ok()) {
      continue;
    }

    const auto sarr = S_new.const_array(mfi);
    const auto parr = pressure_mf->const_array(mfi);
    amrex::ParallelFor(
      bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x =
          prob_lo[0] +
          (static_cast<amrex::Real>(i - ilo) + 0.5) * dx[0];
        const amrex::Real z =
          prob_lo[2] +
          (static_cast<amrex::Real>(k - klo) + 0.5) * dx[2];
        const amrex::Real xr = x - hole_cx;
        const amrex::Real zr = z - hole_cz;
        const bool in_square =
          (xr > -hole_radius) && (xr < hole_radius) &&
          (zr > -hole_radius) && (zr < hole_radius);
        const bool in_circle = xr * xr + zr * zr <= hole_radius2;
        if (
          (inject_fuel == 1 && in_square) ||
          (inject_fuel == 2 && in_circle)) {
          const amrex::Real rho = sarr(i, j, k, URHO);
          const amrex::Real u = sarr(i, j, k, UMX) / rho;
          const amrex::Real v = sarr(i, j, k, UMY) / rho;
          const amrex::Real w = sarr(i, j, k, UMZ) / rho;
          const amrex::Real magvel = std::sqrt(u * u + v * v + w * w);
          const amrex::Real y_h2 = sarr(i, j, k, UFS + H2_ID) / rho;
          const amrex::Real cell_area = dx[0] * dx[2];

          amrex::Gpu::Atomic::Add(&sums[0], parr(i, j, k, 0));
          amrex::Gpu::Atomic::Add(&sums[1], sarr(i, j, k, UTEMP));
          amrex::Gpu::Atomic::Add(&sums[2], u);
          amrex::Gpu::Atomic::Add(&sums[3], v);
          amrex::Gpu::Atomic::Add(&sums[4], w);
          amrex::Gpu::Atomic::Add(&sums[5], magvel);
          amrex::Gpu::Atomic::Add(&sums[6], 1.0);
          amrex::Gpu::Atomic::Add(&sums[7], rho);
          amrex::Gpu::Atomic::Add(&sums[8], cell_area);
          amrex::Gpu::Atomic::Add(&sums[9], rho * v * cell_area);
          amrex::Gpu::Atomic::Add(&sums[10], rho * y_h2 * v * cell_area);
          amrex::Gpu::Atomic::Add(&sums[11], rho * v * v * cell_area);
        }
      });
  }
  amrex::Gpu::streamSynchronize();

  amrex::Vector<amrex::Real> h_sums(nvals, 0.0);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, d_sums.begin(), d_sums.end(), h_sums.begin());
  amrex::ParallelDescriptor::ReduceRealSum(h_sums.data(), h_sums.size());

  if (amrex::ParallelDescriptor::IOProcessor()) {
    const amrex::Real count = h_sums[6];
    const amrex::Real inv_count = (count > 0.0) ? 1.0 / count : 0.0;
    const amrex::Real jet_area = h_sums[8];
    const amrex::Real avg_momentum_flux_y =
      (jet_area > 0.0) ? h_sums[11] / jet_area : 0.0;

    amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
    massfrac_air[O2_ID] = PeleC::h_prob_parm_device->Yox_channel;
    massfrac_air[N2_ID] = 1.0 - PeleC::h_prob_parm_device->Yox_channel;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real rho_air = 0.0;
    amrex::Real cs_air = 0.0;
    eos.PYT2R(
      PeleC::h_prob_parm_device->Pres_amb, massfrac_air,
      PeleC::h_prob_parm_device->Temp_amb, rho_air);
    eos.RTY2Cs(
      rho_air, PeleC::h_prob_parm_device->Temp_amb, massfrac_air, cs_air);

    const amrex::Real u_air =
      PeleC::h_prob_parm_device->Mach_channel * cs_air;
    const amrex::Real crossflow_momentum_flux = rho_air * u_air * u_air;
    const amrex::Real effective_momentum_ratio =
      (crossflow_momentum_flux > 0.0)
        ? avg_momentum_flux_y / crossflow_momentum_flux
        : 0.0;

    const amrex::Real missing =
      std::numeric_limits<amrex::Real>::quiet_NaN();
    write_jet_probe_row(
      step, time, count, jet_area,
      (count > 0.0) ? h_sums[0] * inv_count : missing,
      (count > 0.0) ? h_sums[1] * inv_count : missing,
      (count > 0.0) ? h_sums[7] * inv_count : missing,
      (count > 0.0) ? h_sums[2] * inv_count : missing,
      (count > 0.0) ? h_sums[3] * inv_count : missing,
      (count > 0.0) ? h_sums[4] * inv_count : missing,
      (count > 0.0) ? h_sums[5] * inv_count : missing,
      h_sums[9], h_sums[10], h_sums[11], avg_momentum_flux_y,
      crossflow_momentum_flux, effective_momentum_ratio,
      *PeleC::h_prob_parm_device);
  }
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
}
