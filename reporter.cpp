#include <iostream>

#include "DagMC.hpp"
#include "argparse/argparse.hpp"

int main(int argc, char** argv) {

  argparse::ArgumentParser args("DAGMC Model Reporter");

  args.add_argument("dagmc_filename")
      .help("Path to the DAGMC file.");

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

  return 0;
}
