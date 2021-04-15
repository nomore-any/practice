#pragma once

#include "Element.h"

using namespace std;

//! 从Element类继承的Bar类
class CQ4 : public CElement
{
public:

//!	构造函数
	CQ4();

//!	析构函数
	~CQ4();

//!	读取Q4单元的数据
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	打印Q4单元的数据
	virtual void Write(COutputter& output);

//!	计算单元刚度阵，存到传入的数组里
	virtual void ElementStiffness(double* Matrix);

//!	计算单元应力
	virtual void ElementStress(double* stress, double* Displacement);
};
