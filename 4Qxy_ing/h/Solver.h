/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "SkylineMatrix.h"

//!	������࣬�������KX=F
class CLDLTSolver
{
private:
    
    CSkylineMatrix<double>& K;

public:

//!	���캯��������һ��CSkylineMatrix��Ķ�����г�ʼ��
	CLDLTSolver(CSkylineMatrix<double>* K): K(*K) {};

//!	�Ըն���K����LDLT�ֽ�
	void LDLT();

//!	����ȫ���غ��������ش����λ��
	void BackSubstitution(double* Force); 
};
