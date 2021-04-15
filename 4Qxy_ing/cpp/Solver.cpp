/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Solver.h"

#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>

using namespace std;

// LDLT 分解，分解后的矩阵直接存到原来刚度阵的内存空间(内存复用)。公式见课件2-4第6页
// 最终K矩阵的元素是这样的：对角元存储的是D矩阵的元素，非对角元存储的是L矩阵(上三角阵)的元素
void CLDLTSolver::LDLT()
{
	//从SkylineMatrix类的对象获取全局刚度阵的维数、全局刚度阵的列高
	unsigned int N = K.dim();
    unsigned int* ColumnHeights = K.GetColumnHeights();  

	//因为d11=k11，所以第一列就不用管了，从第二列开始操作
	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n 
	{
        // 计算第j列一个非零元的位置
		unsigned int mj = j - ColumnHeights[j-1];
        
		for (unsigned int i = mj+1; i <= j-1; i++)	
		{
			// 当计算u_ij时，此时K矩阵的元素是什么样的？
			//这个时候j'<j的所有元素都已经是最终的值(L的元素与D的元素)，j'=j,i'<i的元素还是U矩阵的元素


            // 计算第i列第一个非零元位置
			unsigned int mi = i - ColumnHeights[i-1];

			double C = 0.0;
			for (unsigned int r = max(mi, mj); r <= i-1; r++)
				C += K(r,i) * K(r,j);		// C += L_ri * U_rj

			K(i,j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = K(r,j) / K(r,r);	// L_rj = U_rj / D_rr
			K(j,j) -= Lrj * K(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			K(r,j) = Lrj;
		}


		//刚度阵对角元素非正就报错
        if (fabs(K(j,j)) <= FLT_MIN)
        {
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
            	 << "    Euqation no = " << j << endl
            	 << "    Pivot = " << K(j,j) << endl;
            
            exit(4);
        }
    }
};

// 回代求解位移，公式见课件2-4第8页
void CLDLTSolver::BackSubstitution(double* Force)
{
	//从SkylineMatrix对象获取矩阵维数、列高数组
	unsigned int N = K.dim();
    unsigned int* ColumnHeights = K.GetColumnHeights();   

//	Reduce right-hand-side load vector (LV = R)
	for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
        unsigned int mi = i - ColumnHeights[i-1];

		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] -= K(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
	}

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] /= K(i,i);	// Vbar = D^(-1) V

	for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
	{
        unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= K(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
};
