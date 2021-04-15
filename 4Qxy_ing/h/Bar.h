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

#include "Element.h"

using namespace std;

//! 从Element类继承的Bar类
class CBar : public CElement
{
public:

//!	构造函数
	CBar();

//!	析构函数
	~CBar();

//!	读取Bar单元的数据
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	打印Bar单元的数据
	virtual void Write(COutputter& output);

//!	计算单元刚度阵，存到传入的数组里
	virtual void ElementStiffness(double* Matrix);

//!	计算单元应力
	virtual void ElementStress(double* stress, double* Displacement);
};
