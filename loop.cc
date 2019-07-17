#include "connected.h"
#include "loop.h"

// std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> gcoor({24, 24, 24, 64});

using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);

  int t_min = 10; // should be consistent with calculating loop2
  double Mpi = 0.14;

	// int traj_start = 680, traj_num = 70;
  // std::vector<int> trajs(traj_num);
  // for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;
  std::vector<int> trajs {1370};

  cout << "trajs: " << endl;
  cout << trajs << endl;

  for(int traj: trajs) {


    // LatticeLoop loop1(grid);
    // read_loop(loop1, loop_path_24ID(traj));
    //
    // // std::vector<LatticePropagator> wall_props(gcoor[Tdir], grid);
    // // read_wall_src_props("24ID", traj, wall_props);
    // // LatticeLoop loop2 = loops_contraction(wall_props, t_min);
    // // // writeScidac(loop2, "loop2_1370.lat");
    //
    // LatticeLoop loop2(grid);
    // readScidac(loop2, "loop2_1370.lat");
    //
    // LatticePGG disc = three_point_loop(loop1, loop2, t_min, Mpi);

    LatticePGG disc(grid);
    calculate_disc(disc, "24ID", traj, t_min, Mpi);

    cout << disc << endl;
  }


  end();

  return 0;
}
