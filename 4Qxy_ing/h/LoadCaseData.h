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

#include "Outputter.h"

using namespace std;

//! 载荷工况数据类，每增加一个载荷工况就增加一个该类的对象
class CLoadCaseData
{
public:

	unsigned int nloads;	//!该工况的集中载荷数
	unsigned int* node;		//!< 从1~nloads,每个载荷都作用在哪个节点上
	unsigned int* dof;		//!< 从1~nloads,每个载荷都作用在哪个自由度上
	double* load;			//!< 从1~nloads,每个载荷都有多大

public:
	//构造函数，各种东西初始化
	CLoadCaseData() : nloads(0), node(NULL), dof(NULL), load(NULL) {};
	~CLoadCaseData();

//!	传入集中载荷数，在堆区分配空间
	void Allocate(unsigned int num);

//!	读取载荷数据
	bool Read(ifstream& Input);

//!	输出载荷数据
	void Write(COutputter& output);
};
