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

//!	求解器类，用于求解KX=F
class CLDLTSolver
{
private:
    
    CSkylineMatrix<double>& K;

public:

//!	构造函数，传入一个CSkylineMatrix类的对象进行初始化
	CLDLTSolver(CSkylineMatrix<double>* K): K(*K) {};

//!	对刚度阵K进行LDLT分解
	void LDLT();

//!	传入全局载荷向量，回代求解位移
	void BackSubstitution(double* Force); 
};
