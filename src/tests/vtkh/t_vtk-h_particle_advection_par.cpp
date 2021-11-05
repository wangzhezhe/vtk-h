//-----------------------------------------------------------------------------
///
/// file: t_vtk-h_particle_advection_par.cpp
///
//-----------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkh/filters/Streamline.hpp>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/CellSetSingleType.h>
#include "t_test_utils.hpp"
#include <iostream>
#include <mpi.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/filter/PointAverage.h>

void
SaveOutput(vtkh::DataSet* output, int rank, std::string fname)
{
  std::cerr << "rank " << rank << " saving data from " << output->GetNumberOfDomains() << " blocks" << std::endl;
  for (int i = 0; i < output->GetNumberOfDomains(); i++)
  {
    auto ds = output->GetDomain(i);

    char nm[512];
    sprintf(nm, "%s.r%d.%d.vtk", fname.c_str(), rank, i);
    vtkm::io::VTKDataSetWriter wrt(nm);
    wrt.WriteDataSet(ds);
  }
}

void LoadData(std::string& fname, std::vector<vtkm::cont::DataSet>& dataSets, int rank, int nRanks)
{
  std::string buff;
  std::ifstream is;
  is.open(fname);
  std::cout << "Opening: " << fname << std::endl;
  if (!is)
  {
    std::cout << "File not found! : " << fname << std::endl;
    throw "unknown file: " + fname;
  }

  auto p0 = fname.rfind(".visit");
  if (p0 == std::string::npos)
    throw "Only .visit files are supported.";
  auto tmp = fname.substr(0, p0);
  auto p1 = tmp.rfind("/");
  auto dir = tmp.substr(0, p1);

  std::getline(is, buff);
  auto numBlocks = std::stoi(buff.substr(buff.find("!NBLOCKS ") + 9, buff.size()));

  int nPer = numBlocks / nRanks;
  int b0 = rank * nPer, b1 = (rank + 1) * nPer;
  if (rank == (nRanks - 1))
    b1 = numBlocks;

  if (rank == 0)
    std::cout << "numBlocks= " << numBlocks << " "<<b0<<" "<<b1<<std::endl;

  for (int i = 0; i < numBlocks; i++)
  {
    std::getline(is, buff);
    if (i >= b0 && i < b1)
    {
      vtkm::cont::DataSet ds;
      std::string vtkFile = dir + "/" + buff;
      std::cerr << "Rank " << rank << " opening file " << buff << std::endl;
      vtkm::io::VTKDataSetReader reader(vtkFile);
      ds = reader.ReadDataSet();

      vtkm::cont::ArrayHandle<vtkm::FloatDefault> Vx, Vy, Vz;
      Vx = ds.GetField("velocityx").GetData().AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
      Vy = ds.GetField("velocityy").GetData().AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
      Vz = ds.GetField("velocityz").GetData().AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
      auto cellV = vtkm::cont::make_ArrayHandleSOA<vtkm::Vec3f_64>({Vx, Vy, Vz});
      ds.AddField(vtkm::cont::make_FieldCell("velC", cellV));

      vtkm::filter::PointAverage avg;
      avg.SetActiveField("velC");
      avg.SetOutputFieldName("vel");
      ds = avg.Execute(ds);

      //rename ghost field
      auto temp = ds.GetField("topo_ghosts");
      if(temp.GetNumberOfValues() >= 1)
      {
          auto ghostArr = temp.GetData().AsArrayHandle<vtkm::cont::ArrayHandleBasic<vtkm::FloatDefault>>();
          const vtkm::FloatDefault* buff = ghostArr.GetReadPointer();
          vtkm::cont::ArrayHandle<vtkm::UInt8> ghosts;
          ghosts.Allocate(temp.GetNumberOfValues());
          for (vtkm::Id z = 0; z < temp.GetNumberOfValues(); z++)
          {
              ghosts.WritePortal().Set(z, static_cast<vtkm::UInt8>(buff[z]));
          }
          ds.AddCellField("vtkmGhostCells", ghosts);
          //data.GetDomain(i).RemoveField("topo_ghosts");
      }

      /*
      auto f = ds.GetField("grad").GetData();
      vtkm::cont::ArrayHandle<vtkm::Vec<double, 3>> fieldArray;
      f.AsArrayHandle(fieldArray);
      int n = fieldArray.GetNumberOfValues();
      auto portal = fieldArray.WritePortal();
      for (int ii = 0; ii < n; ii++)
        portal.Set(ii, vtkm::Vec<double, 3>(1, 0, 0));
      */

      dataSets.push_back(ds);
    }
  }
}





void checkValidity(vtkh::DataSet *data, const int maxSteps, bool isSL)
{
  int numDomains = data->GetNumberOfDomains();

  //Check all domains
  for(int i = 0; i < numDomains; i++)
  {
    auto currentDomain = data->GetDomain(i);
    auto cs = currentDomain.GetCellSet();
    if (isSL)
    {
      auto cellSet = cs.Cast<vtkm::cont::CellSetExplicit<>>();
      //Ensure that streamlines took <= to the max number of steps
      for(int j = 0; j < cellSet.GetNumberOfCells(); j++)
      {
        EXPECT_LE(cellSet.GetNumberOfPointsInCell(j), maxSteps);
      }
    }
    else
    {
      if (!cs.IsType<vtkm::cont::CellSetSingleType<>>())
        EXPECT_TRUE(false);
    }
  }
}

void writeDataSet(vtkh::DataSet *data, std::string fName, int rank)
{
  int numDomains = data->GetNumberOfDomains();
  std::cerr << "num domains " << numDomains << std::endl;
  for(int i = 0; i < numDomains; i++)
  {
    char fileNm[128];
    sprintf(fileNm, "%s.rank%d.domain%d.vtk", fName.c_str(), rank, i);
    vtkm::io::writer::VTKDataSetWriter write(fileNm);
    write.WriteDataSet(data->GetDomain(i));
  }
}

static inline vtkm::FloatDefault
rand01()
{
  return (vtkm::FloatDefault)rand() / (RAND_MAX+1.0f);
}

static inline vtkm::FloatDefault
randRange(const vtkm::FloatDefault &a, const vtkm::FloatDefault &b)
{
    return a + (b-a)*rand01();
}

template <typename FilterType>
vtkh::DataSet *
RunFilter(vtkh::DataSet& input,
          const std::string& fieldName,
          const std::vector<vtkm::Particle>& seeds,
          int maxAdvSteps,
          double stepSize)
{
  FilterType filter;

  filter.SetInput(&input);
  filter.SetField(fieldName);
  filter.SetNumberOfSteps(maxAdvSteps);
  filter.SetStepSize(stepSize);
  filter.SetSeeds(seeds);
  filter.Update();

  return filter.GetOutput();
}

//----------------------------------------------------------------------------
TEST(vtkh_particle_advection, vtkh_serial_particle_advection)
{
  MPI_Init(NULL, NULL);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  vtkh::SetMPICommHandle(MPI_Comm_c2f(MPI_COMM_WORLD));




 std::string dataFile = "/home/jkress/packages/visualizationPerformanceEvaluation/vtk-h/build/files.visit";
 std::vector<vtkm::cont::DataSet> dataSets;
 LoadData(dataFile, dataSets, rank, size);
 vtkh::DataSet pds;
 for (const auto& ds: dataSets)
   pds.AddDomain(ds, rank);

 vtkh::Streamline pa;
 std::vector<vtkm::Particle> seeds;
 seeds.push_back({{0.0731594, 0.00178438, 0.0588844}, rank});

 if (size == 1)
 {
   std::ofstream out("pts.txt");
   for (const auto& p : seeds)
     out<<p.ID<<", "<<p.Pos[0]<<", "<<p.Pos[1]<<", "<<p.Pos[2]<<std::endl;
 }
 auto seedArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);

 /*
   vtkm::cont::ArrayHandle<vtkm::Particle> seedArray;
   seedArray = vtkm::cont::make_ArrayHandle({ vtkm::Particle(vtkm::Vec3f(.1f, .1f, .9f), 0),
   vtkm::Particle(vtkm::Vec3f(.1f, .6f, .6f), 1),
   vtkm::Particle(vtkm::Vec3f(.1f, .9f, .1f), 2) });
 */
 pa.SetStepSize(0.001f);
 //pa.SetStepSize(0.01f);
 pa.SetNumberOfSteps(1000);
 pa.SetSeeds(seeds);
 pa.SetField("vel");


 pa.SetInput(&pds);
 pa.Update();
 auto output = pa.GetOutput();
 SaveOutput(output, rank, "output");

 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Finalize();



#if 0
  std::cout << "Running parallel Particle Advection, vtkh - with " << comm_size << " ranks" << std::endl;

  vtkh::DataSet data_set;
  const int base_size = 32;
  const int blocks_per_rank = 1;
  const int maxAdvSteps = 1000;
  const int num_blocks = comm_size * blocks_per_rank;

  for(int i = 0; i < blocks_per_rank; ++i)
  {
    int domain_id = rank * blocks_per_rank + i;
    data_set.AddDomain(CreateTestDataRectilinear(domain_id, num_blocks, base_size), domain_id);
  }

  std::vector<vtkm::Particle> seeds;

  vtkm::Bounds bounds = data_set.GetGlobalBounds();
  std::cout<<"Bounds= "<<bounds<<std::endl;

  for (int i = 0; i < 100; i++)
  {
    vtkm::Particle p;
    p.Pos = vtkm::Vec3f(randRange(bounds.X.Min, bounds.X.Max),
                        randRange(bounds.Y.Min, bounds.Y.Max),
                        randRange(bounds.Z.Min, bounds.Z.Max));
    p.ID = static_cast<vtkm::Id>(i);

    seeds.push_back(p);
  }

  vtkh::DataSet *outPA=NULL, *outSL=NULL;

  outPA = RunFilter<vtkh::ParticleAdvection>(data_set, "vector_data_Float64", seeds, maxAdvSteps, 0.1);
  outPA->PrintSummary(std::cout);
  checkValidity(outPA, maxAdvSteps+1, false);

  outSL = RunFilter<vtkh::Streamline>(data_set, "vector_data_Float64", seeds, maxAdvSteps, 0.1);
  outSL->PrintSummary(std::cout);
  checkValidity(outSL, maxAdvSteps+1, true);

  writeDataSet(outSL, "advection_SeedsRandomWhole", rank);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
}
