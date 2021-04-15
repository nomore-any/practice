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

// LDLT �ֽ⣬�ֽ��ľ���ֱ�Ӵ浽ԭ���ն�����ڴ�ռ�(�ڴ渴��)����ʽ���μ�2-4��6ҳ
// ����K�����Ԫ���������ģ��Խ�Ԫ�洢����D�����Ԫ�أ��ǶԽ�Ԫ�洢����L����(��������)��Ԫ��
void CLDLTSolver::LDLT()
{
	//��SkylineMatrix��Ķ����ȡȫ�ָն����ά����ȫ�ָն�����и�
	unsigned int N = K.dim();
    unsigned int* ColumnHeights = K.GetColumnHeights();  

	//��Ϊd11=k11�����Ե�һ�оͲ��ù��ˣ��ӵڶ��п�ʼ����
	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n 
	{
        // �����j��һ������Ԫ��λ��
		unsigned int mj = j - ColumnHeights[j-1];
        
		for (unsigned int i = mj+1; i <= j-1; i++)	
		{
			// ������u_ijʱ����ʱK�����Ԫ����ʲô���ģ�
			//���ʱ��j'<j������Ԫ�ض��Ѿ������յ�ֵ(L��Ԫ����D��Ԫ��)��j'=j,i'<i��Ԫ�ػ���U�����Ԫ��


            // �����i�е�һ������Ԫλ��
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


		//�ն���Խ�Ԫ�ط����ͱ���
        if (fabs(K(j,j)) <= FLT_MIN)
        {
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
            	 << "    Euqation no = " << j << endl
            	 << "    Pivot = " << K(j,j) << endl;
            
            exit(4);
        }
    }
};

// �ش����λ�ƣ���ʽ���μ�2-4��8ҳ
void CLDLTSolver::BackSubstitution(double* Force)
{
	//��SkylineMatrix�����ȡ����ά�����и�����
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
