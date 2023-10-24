#include <iostream>

#include "DagMC.hpp"
#include "argparse/argparse.hpp"


double prn() {
  return (double) std::rand() / (double)RAND_MAX;
}


int main(int argc, char** argv) {

  argparse::ArgumentParser args("DAGMC Model Reporter");

  args.add_argument("dagmc_filename")
      .help("Path to the DAGMC file.");

  args.add_argument("-n")
      .help("Number of find_volume samples to perform")
      .nargs(1)
      .default_value(1000000) // 1M queries by default
      .scan<'i', int>();

  args.add_argument("-s")
      .help("Skip find_volume queries")
      .default_value(false)
      .implicit_value(true);

  try {
    args.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cout << err.what() << std::endl;
    std::cout << args;
    std::exit(1);
  }

  moab::DagMC dag{};
  dag.load_file(args.get<std::string>("dagmc_filename").c_str());
  dag.setup_indices();

  moab::Range tris;
  dag.moab_instance()->get_entities_by_type(0, moab::MBTRI, tris);

  std::cout << "Model contains " << tris.size() << " triangles." << std::endl;

  for (int i = 0; i < dag.num_entities(3); i++) {
    moab::EntityHandle vol = dag.entity_by_index(3, i+1); // Fortran indexing
    double volume;
    dag.measure_volume(vol, volume);
    std::cout << "Volume " << dag.id_by_index(3, i) << ": " << volume << " (cc)" << std::endl;
  }

  if (args.get<bool>("s")) {
    std::cout << "Skipping find_volume queries" << std::endl;
    return 0;
  }

  dag.init_OBBTree();

  // get the bounding box of the entire model
  double lower_left[3], upper_right[3];
  moab::Range vertices;
  dag.moab_instance()->get_entities_by_type(0, moab::MBVERTEX, vertices);

  std::vector<double> x(vertices.size());
  std::vector<double> y(vertices.size());
  std::vector<double> z(vertices.size());
  dag.moab_instance()->get_coords(vertices, x.data(), y.data(), z.data());

  auto xminmax = std::minmax(x.begin(), x.end());
  auto yminmax = std::minmax(y.begin(), y.end());
  auto zminmax = std::minmax(z.begin(), z.end());

  lower_left[0] = *xminmax.first; upper_right[0] = *xminmax.second;
  lower_left[1] = *yminmax.first; upper_right[1] = *yminmax.second;
  lower_left[2] = *zminmax.first; upper_right[2] = *zminmax.second;

  int n_samples = args.get<int>("n");
  moab::EntityHandle result;

  // perform a bunch of find volme queries. This has the effect of randomly sampling ray fires throughout the model
  for (int i = 0; i < n_samples; i++) {
    double xyz[3] = {lower_left[0] + (upper_right[0] - lower_left[0]) * prn(),
                     lower_left[1] + (upper_right[1] - lower_left[1]) * prn(),
                     lower_left[2] + (upper_right[2] - lower_left[2]) * prn()};
    dag.find_volume(xyz, result);
  }


  return 0;
}
