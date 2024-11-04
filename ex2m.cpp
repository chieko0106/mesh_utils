/*                             MFEM modified from ex2
//
// Compile with: make ex2
//
// Sample runs:  ex2 -m ../data/beam-tri.mesh
//               ex2 -m ../data/beam-quad.mesh
//               ex2 -m ../data/beam-tet.mesh
//               ./ex2m -m ../data/beam-hex.mesh
//               ex2 -m ../data/beam-wedge.mesh
//               ex2 -m ../data/beam-quad.mesh -o 3 -sc
//               ex2 -m ../data/beam-quad-nurbs.mesh
//               ./ex2m -m ../data/beam-hex-nurbs.mesh
//
// Description:  This example code solves a simple linear elasticity problem
//               describing a multi-material cantilever beam.
//
//               Specifically, we approximate the weak form of -div(sigma(u))=0
//               where sigma(u)=lambda*div(u)*I+mu*(grad*u+u*grad) is the stress
//               tensor corresponding to displacement field u, and lambda and mu
//               are the material Lame constants. The boundary conditions are
//               u=0 on the fixed part of the boundary with attribute 1, and
//               sigma(u).n=f on the remainder with f being a constant pull down
//               vector on boundary elements with attribute 2, and zero
//               otherwise. The geometry of the domain is assumed to be as
//               follows:
//
//                                 +----------+----------+
//                    boundary --->| material | material |<--- boundary
//                    attribute 1  |    1     |    2     |     attribute 2
//                    (fixed)      +----------+----------+     (pull down)
//
//               The example demonstrates the use of high-order and NURBS vector
//               finite element spaces with the linear elasticity bilinear form,
//               meshes with curved elements, and the definition of piece-wise
//               constant and vector coefficient objects. Static condensation is
//               also illustrated.
//
//               We recommend viewing Example 1 before viewing this example.
*/ 

#include "mfem.hpp"
#include "ex2m.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


int main(int argc, char *argv[])
{
   // 1. 参数设置
   
   const char *mesh_file = "../data/beam-tri.mesh";
   int order = 1;
   int ref_levels = 1; //默认的网格加密层级
   bool static_cond = false;
   bool visualization = false;
   bool paraview=true;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&paraview,"-pa", "--paraview", "-no-pa",
                  "--no-paraview",
                  "if use paraview for visualization.");    
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. 从外部文件读取网格数据
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   if (mesh.attributes.Max() < 2 || mesh.bdr_attributes.Max() < 2)
   {
      cerr << "\nInput mesh should have at least two materials and "
           << "two boundary attributes! (See schematic in ex2.cpp)\n"
           << endl;
      return 3;
   }

   // 3. Select the order of the finite element discretization space. For NURBS
   //    meshes, we increase the order by degree elevation.
   if (mesh.NURBSext)
   {
      mesh.DegreeElevate(order, order);
   }

   // 4. 网格加密
      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }

   // 5. 定义网格空间

   FiniteElementCollection *fec;
   FiniteElementSpace *fespace;
   if (mesh.NURBSext)
   {
      fec = NULL;
      fespace = mesh.GetNodes()->FESpace();
   }
   else
   {
      fec = new H1_FECollection(order, dim);
      fespace = new FiniteElementSpace(&mesh, fec, dim);
   }
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl << "Assembling: " << flush;

   // 6. 设置固定边界条件
   Array<int> ess_tdof_list, ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = 1; // =1 essential 边界条件
   fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   // 7. 设置给定力的边界条件
   VectorArrayCoefficient f(dim);
   for (int i = 0; i < dim-1; i++)
   {
      f.Set(i, new ConstantCoefficient(0.0));
   }
   {
      Vector pull_force(mesh.bdr_attributes.Max()); 
      pull_force = 0.0;
      pull_force(1) = -1.0e-2;// 在bdr标记第2个位置上设置边界力的值，其他边界力为0
      f.Set(dim-1, new PWConstCoefficient(pull_force));//piece-wise constants coefficient
      //f.Set(dim-1,new PWConstCoefficient(pull_force)) 表示物理编号为1的边上沿dim-1(2D-y&3D-z)方向施加一个大小为-1.0e-2的力。
   }

   LinearForm *b = new LinearForm(fespace);
   b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
   cout << "r.h.s. ... " << flush;
   b->Assemble();

   // 8. 设置初始条件
   GridFunction x(fespace);
   x = 0.0;

   // 9. 对不同的区域设定不同的材料参数
   Vector lambda(mesh.attributes.Max());
   lambda = 1.0;
   lambda(0) = lambda(1)*50;
      //向量lambda = {50 1}表示实体区域0上拉梅常数的lamda为50而实体区域1则为1
   PWConstCoefficient lambda_func(lambda);//
      //PWConstCoefficient lambda_func(lambda) 表示实体区域0上所有单元定义值为50的拉梅常数而实体区域1所有单元值为1。
   Vector mu(mesh.attributes.Max());
   mu = 1.0;
   mu(0) = mu(1)*50;
   PWConstCoefficient mu_func(mu);
   BilinearForm *a = new BilinearForm(fespace);//定义刚度矩阵
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));

   // 10. 集成刚度矩阵
   cout << "matrix ... " << flush;
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();//生成刚度矩阵
   SparseMatrix A;
   Vector B, X;
   //消元
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   //小写的a，x,*b指消元前的刚度矩阵、未知量和右端项， A, X, B 指消元后刚度矩阵、未知量和右端项。同时需要提供根据标记向量ess_bdr确定的边界节点自由度标记数组ess_tdof_list
   cout << "done." << endl;

   cout << "Size of linear system: " << A.Height() << endl;

#ifndef MFEM_USE_SUITESPARSE
   // 11. PCG求解器
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 500, 1e-8, 0.0);
#else
   // 11. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(B, X);
#endif

   // 12. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, x);

   // 13. For non-NURBS meshes
   if (!mesh.NURBSext)
   {
      mesh.SetNodalFESpace(fespace);
   }
   // 计算应变
   int vdim = dim*(dim+1)/2;
   FiniteElementSpace *fespaceE = new FiniteElementSpace(&mesh, fec,vdim ); 
   GridFunction strain_gf(fespaceE);   
   StrainCoefficient *strain_coef = new StrainCoefficient(vdim); 
   strain_coef->SetDisplacement(x);
   strain_gf.ProjectCoefficient(*strain_coef);  //仅调用ProjectCoefficient函数即完成投影，这个过程引入投影误差
   // 计算应力
   GridFunction stress_gf(fespaceE);  
   StressCoefficient *stress_coef = new StressCoefficient(vdim,lambda_func, mu_func);
   stress_coef->SetDisplacement(x);
   stress_gf.ProjectCoefficient(*stress_coef);



   // 15. 保存可视化结果
   mfem::ParaViewDataCollection paraview_dc("ex2m", &mesh);
   if (paraview)
   {
      paraview_dc.SetPrefixPath("ParaView");
      paraview_dc.SetLevelsOfDetail(order);
      paraview_dc.SetDataFormat(VTKFormat::BINARY);
      paraview_dc.SetHighOrderOutput(true);
      paraview_dc.SetCycle(0);
      paraview_dc.SetTime(0.0);
      paraview_dc.RegisterField("displacement",&x);
      //paraview_dc.RegisterField("stress_normal",&stress_normal);
      //paraview_dc.RegisterField("stress_shear",&stress_shear);
      paraview_dc.RegisterField("strain",&strain_gf);
      paraview_dc.RegisterField("stress",&stress_gf);
      //paraview_dc.RegisterField("Pressure",&p_gf);
      paraview_dc.Save();
   }
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << mesh << x << flush;
   }

   // 16. Free the used memory.
   delete a;
   delete b;
   if (fec)
   {
      delete fespace;
      delete fec;
   }

   return 0;
}

