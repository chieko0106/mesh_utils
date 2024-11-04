#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace mfem;

// 后处理-应变
class StrainCoefficient : public VectorCoefficient    //定义向量系数的子类，明确继承关系。
{
protected:
   GridFunction *u; // 要传入的位移场
   DenseMatrix gshape; //单元内所有节点的形函数的梯度组成的矩阵 
public:
   StrainCoefficient(int vdim_)
    : VectorCoefficient(vdim_) { }  //构造函数，指定VectorCoefficient的分量数，vdim = dim*(dim+1)/2 。
   void SetDisplacement(GridFunction &u_) { u = &u_;}  
   virtual void Eval(Vector &, 
                     ElementTransformation &, 
                     const IntegrationPoint &);  //定义积分评估函数
};

// 后处理-米塞斯应力

class StressCoefficient : public VectorCoefficient
{
protected:
   Coefficient &lambda, &mu;
   GridFunction *u; // displacement
   DenseMatrix grad; // auxiliary matrix, used in Eval 
public:
   StressCoefficient (int vdim_, Coefficient &lambda_, Coefficient &mu_)
    : VectorCoefficient(vdim_),lambda(lambda_),mu(mu_), u(NULL) { }
   void SetDisplacement(GridFunction &u_) { u = &u_;}  
   virtual void Eval(Vector &, ElementTransformation &, const IntegrationPoint &);
};


void StrainCoefficient::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "displacement field is not set");  //如果主程序中
                                   //未向类StrainCoefficient 中传入位移场，则报错并终止主程序。
   const FiniteElementSpace *fes=u->FESpace ();          //根据位移场确定有限元空间
   const FiniteElement *FElem = fes->GetFE(T.ElementNo); //根据形函数类确定当前单元，
                                                         //FElem用于确定单元有关的信息，比如节点数。
   int dim = fes->GetMesh()->Dimension();
   int vdim = dim*(dim+1)/2;       //3D-vdim=6&2D-vdim=3
   int dof = FElem->GetDof();      //dof-单元内节点数。
   Vector displacement(dim*dof);   //单元内位移场节点值的向量。
   Array<int> dofs;                //自由度标记数组申明。
   fes->GetElementDofs(T.ElementNo, dofs);   
   fes->DofsToVDofs(dofs);  //根据当前单元的序数T.ElementNo确定该单元各节点的自由度序列，存入dofs中（默认所有分量）。
   u->GetSubVector(dofs, displacement);      //将单元内位移场的节点值从所有节点位移值u分离出来存入displacement。
                                   //displacement = {u1_1,...,u1_dof,u2_1,...,u2_dof,u3_1,...,u3_dof}注意分量排布方式。
   V.SetSize(vdim);                //输出的积分点应变存入V。
   gshape.SetSize(dof,dim);        //存储节点基函数梯度在积分点上的值的矩阵，是一个dof*dim的矩阵。
   T.SetIntPoint(&ip);             //将积分点ip代入单元所有形函数中。
   FElem->CalcPhysDShape(T,gshape);//得到存储节点基函数梯度在积分点上的值的矩阵gshape。  
   
   V = 0.0;
   if(dim == 3)   //3D strain = {e11,e22,e33,e12,e13,e23}
   {
      for(int i=0; i<dof; i++)
      {
          V(0) += gshape(i,0)*displacement(i); //e11 = n,1_k * u1_k   注意角标爱因斯坦求和。
          V(1) += gshape(i,1)*displacement(i+dof); //e22 = n,2_k * u2_k
          V(2) += gshape(i,2)*displacement(i+dof*2); //e33 = n,3_k * u3_k
          V(3) += 0.5 * (gshape(i,1)*displacement(i) + gshape(i,0)*displacement(i+dof)); 
          //e12 = 0.5*(u1_k * n,2_k + u2_k * n,1_k)
          V(4) += 0.5 * (gshape(i,2)*displacement(i) + gshape(i,0)*displacement(i+dof*2)); 
          //e13 = 0.5*(u1_k * n,3_k + u3_k * n,1_k)
          V(5) += 0.5 * (gshape(i,1)*displacement(i+dof*2) + gshape(i,2)*displacement(i+dof)); 
          //e23 = 0.5*(u3_k * n,2_k + u2_k * n,3_k)
      }
    }else if(dim == 2){ //2D strain = {e11,e22,e12}
       for(int i=0; i<dof; i++)
       {
          V(0) += gshape(i,0)*displacement(i);     //e11 = n,1_k * u1_k  
          V(1) += gshape(i,1)*displacement(i+dof); //e22 = n,2_k * u2_k
          V(2) += 0.5 * (gshape(i,1)*displacement(i) + gshape(i,0)*displacement(i+dof)); 
          //e12 = 0.5*(u1_k * n,2_k + u2_k * n,1_k)
       }
    }  
}


void StressCoefficient::Eval(Vector &Stress,ElementTransformation &T, const IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "displacement field is not set");
   double L = lambda.Eval(T, ip);
   double M = mu.Eval(T, ip);
   int dim=u->FESpace()->GetMesh()->Dimension();
   int vdim = (dim+1)*dim/2;
   T.SetIntPoint(&ip);
   u->GetVectorGradient(T, grad);   
   double div_u = grad.Trace();
   //Vector Stress;
   Stress.SetSize(vdim);
   for(int i=0;i<dim;i++)
   	   for(int j=i;j<dim;j++)
   	   {
   	   	 if (i == j)
   	   	 	Stress(i)=L*div_u + 2*M*grad(i,i);
   	   	 else 
   	   	 	Stress(i+j+dim-1)=M*(grad(i,j) + grad(j,i));  
   	   } 
   
  double MisesStress;  
  if(vdim==6){
      MisesStress = sqrt(0.5*( pow(Stress(0)-Stress(1),2) + pow(Stress(1)-Stress(2),2) 
                  + pow(Stress(0)-Stress(2),2) + 6*( pow(Stress(3),2) + pow(Stress(4),2) + pow(Stress(5),2))));
  }
  else
  {// mises stress 2D = sqrt(s11^2 + s22^2 - s11s22 + 3*s12^2)
      MisesStress = sqrt(pow(Stress(0),2)+pow(Stress(1),2)-Stress(0)*Stress(1)+3*pow(Stress(2),2));
  } 

}

