#include "prob.H"

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

// Read a PMF / Tecplot POINT file (X, temp, u, rho, species mole fractions...)
// and stash the data on host + device, mirroring Exec/RegTests/PMF.
void
read_pmf(const std::string& myfile)
{
  std::string firstline;
  std::string secondline;
  std::string remaininglines;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  PeleC::h_prob_parm_device->pmf_N = line_count;
  PeleC::h_prob_parm_device->pmf_M = variable_count - 1;
  PeleC::prob_parm_host->h_pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->h_pmf_Y.resize(
    static_cast<long>(PeleC::h_prob_parm_device->pmf_N) *
    PeleC::h_prob_parm_device->pmf_M);
  PeleC::prob_parm_host->pmf_Y.resize(
    static_cast<long>(PeleC::h_prob_parm_device->pmf_N) *
    PeleC::h_prob_parm_device->pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (int i = 0; i < PeleC::h_prob_parm_device->pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->h_pmf_X[i];
    for (int j = 0; j < PeleC::h_prob_parm_device->pmf_M; j++) {
      sinput >> PeleC::prob_parm_host
                  ->h_pmf_Y[j * PeleC::h_prob_parm_device->pmf_N + i];
    }
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_X.begin(),
    PeleC::prob_parm_host->h_pmf_X.end(), PeleC::prob_parm_host->pmf_X.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_Y.begin(),
    PeleC::prob_parm_host->h_pmf_Y.end(), PeleC::prob_parm_host->pmf_Y.begin());
  PeleC::h_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::h_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
}

// Locate the flame front in the reference profile (max |dT/dx|, where the
// temperature column is variable j = 0, i.e. h_pmf_Y[i]).
amrex::Real
find_reference_front()
{
  const int N = PeleC::h_prob_parm_device->pmf_N;
  const auto& X = PeleC::prob_parm_host->h_pmf_X;
  const auto& Y = PeleC::prob_parm_host->h_pmf_Y; // temp is the first variable
  amrex::Real best_grad = -1.0;
  amrex::Real x_front = X[N / 2];
  for (int i = 0; i < N - 1; i++) {
    const amrex::Real g = std::abs((Y[i + 1] - Y[i]) / (X[i + 1] - X[i]));
    if (g > best_grad) {
      best_grad = g;
      x_front = 0.5 * (X[i] + X[i + 1]);
    }
  }
  return x_front;
}

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
  std::string pmf_datafile = "LiDryer_H2_p1_phi0_4000tu0300.dat";

  amrex::ParmParse pp("prob");
  pp.query("pamb", PeleC::h_prob_parm_device->pamb);
  pp.query("x_flame", PeleC::h_prob_parm_device->x_flame);
  pp.query("pmf_datafile", pmf_datafile);
  pp.query("u_adv", PeleC::h_prob_parm_device->u_adv);
  if (PeleC::h_prob_parm_device->u_adv >= 0.0) {
    amrex::Print() << "CLEM1DMach: HIGH-SPEED ADVECTION - carrier velocity "
                   << "forced to " << PeleC::h_prob_parm_device->u_adv
                   << " cm/s (prob.u_adv)" << std::endl;
  }

  read_pmf(pmf_datafile);

  // Auto-place the flame: query the profile at (x - standoff), so the
  // reference front (file coordinate x_ref) appears at x = x_flame.
  const amrex::Real x_ref = find_reference_front();
  PeleC::h_prob_parm_device->standoff =
    PeleC::h_prob_parm_device->x_flame - x_ref;
  pp.query("standoff", PeleC::h_prob_parm_device->standoff);

  amrex::Print() << "CLEM1DFlame: reference front at x_ref = " << x_ref
                 << " cm, placed at x_flame = "
                 << PeleC::h_prob_parm_device->x_flame
                 << " cm (standoff = " << PeleC::h_prob_parm_device->standoff
                 << ")" << std::endl;
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
